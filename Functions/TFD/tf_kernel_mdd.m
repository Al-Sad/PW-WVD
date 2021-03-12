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
%
% This is a modified implementation of the original code found in the
% supplementary material of Boashash, Boualem, and Samir Ouelha. "Designing
% high-resolution time–frequency and time–scale distributions for the
% analysis and classification of non-stationary signals: a tutorial review
% with a comparison of features performance." Digital Signal Processing 77
% (2018): 120-152.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                The multi-directional Doppler-lag kernel
%
%  Syntax : G = tf_kernel_mdd(N,af,c,thr)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% N   : Total number of time samples.
% af  : The signal Ambiguity function.
% c   : The kernel shape parameter.
% thr : The kernel threshold from 0 to 1.
%
% <OUTPUTs>
% G : The multi-directional Doppler-lag kernel (N x N).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function G = tf_kernel_mdd(N,af,c,thr)
%% Parameters
LL = 10;
K  = 91;
D0 = 0.01;

%% Computing LFM angles using Radon transform
af = af';
T = linspace(-90,90,K);
[R,xp] = radon(abs(af./max(max(abs(af)))),T);
R  = R(xp==0,:);
Id = speye(K);
D2 = spdiags(ones(K-2,1)*[1 -2 1],0:2, K-2,K);
R  = (Id-inv(Id+LL^2*(D2'*D2)))*R';
R  = R/max(R);

%% Detecting LFM angles
sa = sign(diff([-inf ; R]));
sb = sign(diff([-inf ; R(end:-1:1)]));
sb = sb(end:-1:1);
I  = find(sa==1 & sb==1);
P  = R(sa==1 & sb==1);
T  = 90 - T(I(P > thr));
E  = P(P > thr)/2;
D  = D0*ones(size(T));

%% MDD Kernel in the Doppler-Lag Domain
T = T*pi/180;
y = (-0.5+1/N:1/N:0.5-1/N);
x = (-0.5+1/N:1/N:0.5-1/N);
G = zeros(N,N);
for t = 1:length(T)
    for i = 1:N-1
        for j = 1:N-1
            x1 = cos(T(t))*x(j) - sin(T(t))*y(i);
            y1 = sin(T(t))*x(j) + cos(T(t))*y(i);           
            if(abs(x1) < D(t) && abs(y1) < E(t))
                G(j,i) = G(j,i) + (exp(c*(1-exp(abs(x1/D(t)).^2)))).*(exp(c*(1-exp(abs(y1/E(t)).^2))));
            end
        end
    end
end
G = G'/max(max(G));
end