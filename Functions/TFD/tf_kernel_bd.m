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
%                  The B-distribution Doppler-lag kernel
%
%  Syntax : G = tf_kernel_bd(N,alpha)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% N     : Total number of time samples.
% alpha : The kernel smoothing parameter ranging from 0 to 1.
%
% <OUTPUTs>
% G : The B-distribution Doppler-lag kernel (N x N).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function G = tf_kernel_bd(N,alpha)
%% Generate the BD Doppler-Lag filter G(\tau,\upsilon)
t   = linspace(-N/2,N/2,N);
G1  = cosh(t).^(-2*alpha);
G1  = G1./sum(G1);
G1  = abs(fftshift(fft(G1)));
G1  = repmat(G1,N,1);
Lh  = fix((N-1)/2);
tau = linspace(-Lh,Lh,N);
g2  = abs(tau).^alpha;
g2  = repmat(g2,N,1);
G   = G1'.*g2;
G   = G./max(max(G));
G   = fftshift(G,2);
end