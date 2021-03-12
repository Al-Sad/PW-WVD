% Copyright (c) 2021 Mohammad Fathi Al-Sa'd
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% Email: mohammad.al-sad@tuni.fi, alsad.mohamed@gmail.com
%
% The following reference should be cited whenever this script is used:
% M. Al-Sa'd, B. Boashash, M. Gabbouj, "Design of an Optimal Piece-Wise
% Spline Wigner-Ville Distribution for TFD Performance Evaluation and
% Comparison", IEEE Transactions on Signal Processing, 2021.
%
% Last Modification: 12-March-2021
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%             Mono-Component Non-Stationary Signal Parameters
%
%  Syntax : [IF,IP,IA,w,ptf] = signal_parameters(R,fs,N,M)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% R  : The IF ploynomial order.
% fs : Sampling frequency in Hz.
% N  : Total number of time samples.
% M  : Total number of frequency samples.
%
% <OUTPUTs>
% IF  : The signal IF (1 x N).
% IP  : The signal IP (1 x N).
% IA  : The signal IA (1 x N).
% w   : The signal time-support (1 x N).
% ptf : The IF turning points.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [IF,IP,IA,w,ptf] = signal_parameters(R,fs,N,M)
%% Parameters
t    = 0:1/fs:(N-1)/fs;    % Time array in seocnds
f    = 0:fs/(2*M-1):fs/2;  % Frequency array in hertz
Lmin = 0.1*sqrt(N*M);      % IF minimum Euclidean length
Lmax = 0.9*sqrt(N*M);      % IF maximum Euclidean length
Fmin = 0.1*fs/2;           % IF minimum frequency
Fmax = 0.9*fs/2;           % IF maximum frequency

%% Generate time arrays to the power 0:R and 1:(R+1)
tif = zeros(N,R+1);
tif(:,end) = ones(N,1);
for r = 1:R
    tif(:,r) = t.^(R-r+1);
end
tip = [(t.^(R+1))' tif(:,1:R)]./((R+1):-1:1);
RR = repmat(R:-1:0,R+1,1);

%% Generate polynomial IF of order R alongside its IP
repeat_cycle = 1;
while(repeat_cycle)
    % Generate non-repeating random time points
    pt = t(randperm(N,R+1))';
    w = (t >= min(pt) & t <= max(pt));
    % Generate random frequency points between Fmin and Fmax
    pf = f(randi([round(0.1*M) round(0.9*M)],[R+1 1]))';
    % Create the system matrix A
    A = pt.^RR;
    % Check the matrix rank to avoid dependicies and errors
    if(rank(A) < (R+1))
        repeat_cycle = 1;
    else
        % Create the R order polynomial IF within the lattice T x fs/2
        C  = A\pf;
        IF = (tif*C)'.*w;
        IF(1,~w) = nan;
        % Validate the polynomial order using polynomial fitting
        if(R>1)
            [~, gof] = fit(t(w)', IF(w)', fittype(['poly' num2str(R-1)]),'Normalize','on');
            r2 = gof.adjrsquare;
        else
            r2 = 0;
        end
        % Create the R+1 order polynomial IP
        if(r2 >= 0.95)
            repeat_cycle = 1;
        else
            IP = (tip*C)'.*w;
            IP(1,~w) = 0;
            % Check boundaries and IF curve minimum length
            if((max(IF) <= Fmax && min(IF) >= Fmin))
                nmin = round(min(pt)*fs+1);
                nmax = round(max(pt)*fs+1);
                L = sum(sqrt(1 + diff(IF(nmin:nmax)).^2));
                repeat_cycle = (L < Lmin) || (L > Lmax);
            end
        end
    end
end
ptf = [pt pf];

%% Generate random IA between 0.5 and 1
a = 0.5*rand(1) + 0.5;
IA = repelem(a,1,N);
IA(1,~w) = 0;
end