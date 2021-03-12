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
%             Multi-Component Non-Stationary Signal Generator
%
%  Syntax : [z,IF,IP,IA,Wp,ptf] = signal_generator(R,fs,N,M)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% R  : The IF ploynomial order (1 x P) where P is the number of components.
% fs : Sampling frequency in Hz.
% N  : Total number of time samples.
% M  : Total number of frequency samples.
%
% <OUTPUTs>
% z   : The generated signal (1 x N).
% IF  : The signal IF (P x N).
% IP  : The signal IP (P x N).
% IA  : The signal IA (P x N).
% Wp  : The signal time-support (P x N).
% ptf : The IF turning points.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [z,IF,IP,IA,Wp,ptf] = signal_generator(R,fs,N,M)
%% Initialization
z   = zeros(1,N);
P   = length(R);
IF  = zeros(P,N);
IP  = zeros(P,N);
IA  = zeros(P,N);
Wp  = zeros(P,N);
ptf = cell(P,1);

%% Generate multi-component non-stationary signal z(t)
for p = 1:P
    [IF(p,:), IP(p,:), IA(p,:), Wp(p,:), ptf{p}] = signal_parameters(R(p),fs,N,M);   
    z = z + IA(p,:).*exp(1j*2*pi*IP(p,:)).*Wp(p,:);
end
end
