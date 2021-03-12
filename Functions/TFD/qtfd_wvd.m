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
%                      The Wigner-Ville distribution
%
%  Syntax : TF = qtfd_wvd(s,L)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% s : The input signal (1 x N).
% L : The lag window odd length.
%
% <OUTPUTs>
% TF : The Wigner-Ville distribution of the input signal (N x N).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function TF = qtfd_wvd(s,L)
%% Initialization
N = length(s);
Lh  = fix(L/2);
tau = -Lh:Lh;
mm  = 1 + rem(N + tau, N);
TF  = zeros(N,N);

%% Computing the Analytic Associate of the Input Signal s
z = fft(s,N);           % Transforming S1 into the frequency domain
z(1) = z(1);            % Holding the zero frequency amplitude/energy
z(2:N/2) = 2*z(2:N/2);  % Accepting and doubling half of the frequencies
z(N/2+1:N) = 0;         % Rejecting half of the frequencies
z = ifft(z);            % Transforming Z into the time domain
z = [z zeros(1,N)];     % Padding with zeros

%% Computing the Signal Kernel in the Time-Lag Domain
K = 2*N;
K_TL = zeros(K, K);
for n = 1:N
    Z1 = z(1 + rem(K + n - 1 + tau, K));
    Z2 = z(1 + rem(K + n - 1 - tau, K));
    K_TL(n,mm) = Z1.*conj(Z2);
end

%% Computing the Wigner-Ville Distribution
temp = fft(K_TL,N,2);
TF(1:N-1,:) = real(temp(1:N-1,:));
end