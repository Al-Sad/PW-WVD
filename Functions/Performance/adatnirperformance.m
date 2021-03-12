function [P,AF,AT] = adatnirperformance(W,Tstart,Tend,tol)
% Automatically detects auto-terms in time-frequency representations of
% two-component signals according to Reinhold's ADAT and calculates the 
% normalised instantaneous resolution performance measure.
%
% Output
% P         Numerical performance measure
% AF        Frequency locations of auto-terms, in samples, matrix with two
%           columns
% AT        Time locations of auto-terms, in samples, matrix with two
%           columns
%
% Input
% W         Time-frequency representation, FxT real-valued matrix
% Tstart    First time slice used to calculate performance, default = 1
% Tend      Last time slice used to calculate performance, default = T
% tol       Small error tolerance in estimation of frequency separation,
%           in samples, default = 3
%
% Implemented by Isabella Reinhold, 2016

if nargin < 2 || isempty(Tstart)
    Tstart = 1;
end
if nargin < 3 || isempty(Tend)
    Tend = size(W,2);
end
if nargin < 4 || isempty(tol)
    tol = 3;
end

P = 0;
AF = zeros(Tend-Tstart+1,2);
AT = repmat((Tstart:Tend)',1,2);
k = 1;
for t0 = Tstart:Tend
    if(sum(W(:,t0))>0)
        [Pt0,AF(k,1),AF(k,2)] = reinholdadat(W(:,t0),tol);
        P = P + Pt0;
        k = k + 1;
    end
end
P = P / (Tend - Tstart + 1);
end


function [Pt0,aloc1,aloc2] = reinholdadat(Wt0,tol)

Wt0 = Wt0 / max(abs(Wt0));

% Lamda
[pks,~] = findpeaks(Wt0,'SORTSTR','descend','NPEAKS',2);
lambda = min(pks)/2;

% Delta
[~,locs_all] = findpeaks(Wt0,'MINPEAKHEIGHT',lambda);
b1 = max(1,locs_all(1) - tol);
b2 = min(length(Wt0),locs_all(end) + tol);
delta = max(2*tol,round((b2-b1-2*tol)/2-2*tol));

% Auto-terms
poorRes = 0;
[pks,locs] = findpeaks(Wt0(b1:b2),'SORTSTR','descend','MINPEAKHEIGHT',lambda,'MINPEAKDISTANCE',delta);
if length(locs) > 1
    locs = locs + b1 - 1;
    [locs,I] = sort(locs,'ascend');
    pks = pks(I);
    aloc1 = locs(1);
    apk1 = pks(1);
    aloc2 = locs(end);
    apk2 = pks(end);
else
    poorRes = 1;
    aloc1 = 0;
    aloc2 = 0;
end

if poorRes
    Pt0 = 0;
else
    % Cross-term
    mid = round((aloc1+aloc2)/2);
    tol2 = round((aloc2-aloc1)/8);
    [cpk_pos,cloc_pos] = findpeaks(Wt0(mid-tol2:mid+tol2),'SORTSTR','descend','NPEAKS',1);
    [cpk_neg,cloc_neg] = findpeaks(-Wt0(mid-tol2:mid+tol2),'SORTSTR','descend','NPEAKS',1);
    if ~isempty(cpk_neg) && ~isempty(cpk_pos)
        if cpk_pos > cpk_neg
            cpk = cpk_pos;
            cloc = cloc_pos + mid - tol2 - 1;
        else
            cpk = cpk_neg;
            cloc = cloc_neg + mid - tol2 - 1;
        end
    elseif ~isempty(cpk_pos)
        cpk = cpk_pos;
        cloc = cloc_pos + mid - tol2 - 1;
    elseif ~isempty(cpk_neg)
        cpk = cpk_neg;
        cloc = cloc_neg + mid - tol2 - 1;
    else
        cloc = mid;
        cpk = abs(Wt0(cloc));
    end
    
    % Sidelobes
    if aloc1 >= 3
        spk1Left = findpeaks(abs(Wt0(1:aloc1)),'SORTSTR','descend','NPEAKS',1);
    else
        spk1Left = 0;
    end
    spk1_right = findpeaks(abs(Wt0(aloc1:cloc)),'SORTSTR','descend','NPEAKS',1);
    if ~isempty(spk1_right) && ~isempty(spk1Left)
        spk1 = max(spk1_right,spk1Left);
    elseif ~isempty(spk1_right)
        spk1 = spk1_right;
    elseif ~isempty(spk1Left)
        spk1 = spk1Left;
    else
        spk1 = 0;
    end
    
    if length(Wt0) - aloc2 >= 2
        spk2Right = findpeaks(abs(Wt0(aloc2:end)),'SORTSTR','descend','NPEAKS',1);
    else
        spk2Right = 0;
    end
    spk2_left = findpeaks(abs(Wt0(cloc:aloc2)),'SORTSTR','descend','NPEAKS',1);
    if ~isempty(spk2_left) && ~isempty(spk2Right)
        spk2 = max(spk2_left,spk2Right);
    elseif ~isempty(spk2_left)
        spk2 = spk2_left;
    elseif ~isempty(spk2Right)
        spk2 = spk2Right;
    else
        spk2 = 0;
    end
    
    % Performance
    Pt0 = nirt0(Wt0,[apk1 apk2],[aloc1 aloc2],cpk,[spk1 spk2]); 
end
end

function Pt0 = nirt0(Wt0,apk,aloc,cpk,spk)

Am = (apk(1) + apk(2))/2;
As = (spk(1) + spk(2))/2;
Ax = cpk;

max_diff = round((aloc(2) - aloc(1))/4);

start = max(1,aloc(1)-max_diff);
fin = min(length(Wt0),aloc(1)+max_diff);
W_interest = Wt0(start:fin);
V1 = findv(W_interest,apk(1)*sqrt(2)/2);

start = max(1,aloc(2)-max_diff);
fin = min(length(Wt0),aloc(2)+max_diff);
W_interest = Wt0(start:fin);
V2 = findv(W_interest,apk(2)*sqrt(2)/2);
V = (V1 + V2)/2;

D = 1 - V/(aloc(2) - aloc(1));

Pt0 = 1 - 1/(3)*(As/Am + 1/2*Ax/Am + (1-D));
end

function V = findv(Wt0,mag)

temp = inf;
for i = 1:length(Wt0)
    if abs(Wt0(i) - mag) < temp
        temp = abs(Wt0(i) - mag);
        n = i;
    end
end

a = find(Wt0 == max(Wt0));
V = abs(a - n)*2;
end
