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
%                              The Ideal TFD
%
%  Syntax : [Iz, Ip] = ideal_tfd(IF,IA,W)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% IF : The signal IF (P x N).
% IA : The signal IA (P x N).
% W  : The signal time-support (P x N).
%
% <OUTPUTs>
% Iz : Ideal TFD of the signal (N x N).
% Ip : Ideal TFD of each signal component.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [Iz, Ip] = ideal_tfd(IF,IA,W)
%% Initialization
[P,N] = size(IF);
W     = logical(W);
Lh    = fix((N-1)/2);
tau   = -Lh:Lh;
mm    = 1 + rem(N + tau, N);
K     = 2*N;
n     = 1:N;
K_TL  = zeros(K, K);
Ip    = cell(1,P);

%% Main
Iz = zeros(N,N);
for p = 1:P
    IF(p,~W(p,:)) = nan; % Remove the IF samples that are outside the time support
    wTF = zeros(N,N);
    TFp = zeros(N,N);
    z = [IA(p,:).*W(p,:) zeros(1,N)];
    
    for i = 1:N
        Z1 = z(1 + rem(K + i - 1 + tau, K));
        Z2 = z(1 + rem(K + i - 1 - tau, K));
        K_TL(i,mm) = Z1.*conj(Z2);
    end
    temp = real(fft(K_TL,N,2));
    wTF(1:N-1,:) = fftshift(temp(1:N-1,:),2);
    
    nn = n(W(p,:));
    for i = 1:length(nn)
        IFdiff = IF(p,nn(i)) - 0.25;
        Ndiff  = round(2*N*IFdiff);
        if(Ndiff>=0)
            TFp(nn(i),:) = [zeros(1,Ndiff) wTF(nn(i),1:end-Ndiff)];
        else
            TFp(nn(i),:) = [wTF(nn(i),(abs(Ndiff)+1):end) zeros(1,abs(Ndiff))];
        end
    end
    Ip{p} = TFp;
    Iz = Iz + Ip{p};
end
end