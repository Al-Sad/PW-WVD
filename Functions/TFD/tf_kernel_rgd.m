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
%                 The radial Gaussian Doppler-lag kernel
%
%  Syntax : G = tf_kernel_rgd(N,af,alpha)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% N     : Total number of time samples.
% af    : The signal Ambiguity function.
% alpha : The kernel smoothing parameter.
%
% <OUTPUTs>
% G : The radial Gaussian Doppler-lag kernel (N x N).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function G = tf_kernel_rgd(N,af,alpha)
%% Parameters
muk      = 1; 			   % Step size - controls how fast we climb
eps1     = 0.1;            % Step size test parameter 1
eps2     = 1-eps1;         % Step size test parameter 2
delta    = 10;             % Step size mulitplicative inc/dec-rement
M        = 1001;           % Largest possible step size
eps3     = 1e-5; 		   % Controls how soon we quit climbing
max_iter = 50 ;            % Number of iterations before we stop

%% Convert AF to polar coordinates
tau = linspace(-sqrt(pi*N),(N-1)*sqrt(pi/N),N);
TAU = tau(ones(N,1),:)';
theta = linspace(-sqrt(pi*N),(N-1)*sqrt(pi/N),N);
THETA = theta(ones(N,1),:);
P = round(N/sqrt(2));
rho = linspace(0,sqrt(2*pi*N),P);
dr = 2*sqrt(pi/N);
Q = N;
psi = linspace(0,pi,Q);
dpsi = psi(2)-psi(1);
X = rho'*cos(psi);
Y = rho'*sin(psi);
polamb2 = interp2(THETA,TAU,abs(af),X,Y);
outsquare = find(isnan(polamb2));
polamb2(outsquare) = zeros(size(outsquare));
polamb2 = abs(polamb2/max(max(abs(polamb2)))).^2;

%% RGD Kernel in the Doppler-Lag Domain
gamma  = sqrt(2*pi*alpha/dpsi) ;
sigma1 = ones(1,Q)*(gamma/sqrt(Q)) ;
sigma2 = sigma1 ;
Pindices = 0:P-1 ;
P1  = (Pindices(ones(1,Q),:)).' ;
P3  = (Pindices(ones(1,Q),2:P).^3).' ;
criterion1 = inf ;
criterion2 = inf ;
n_iter = 0 ;
Phi2 = exp(-((rho').^2)*(sigma2.^2).^(-1)) ;
f2 = sum(sum(P1.*polamb2.*Phi2)) ;
while ((criterion1>0) || (criterion2>0)) && n_iter<=max_iter
    sigma1 = sigma2 ;
    f1 = f2 ;
    Phi1 = exp(-((rho').^2)*(sigma1.^2).^(-1)) ;
    grad = sum(P3.*polamb2(2:P,:).*Phi1(2:P,:)).*(2*dr^2*(sigma1).^(-3)) ;
    sigma2 = (sigma1 + muk*grad).*gamma./norm(sigma1 + muk*grad) ;
    Phi2 = exp(-((rho').^2)*(sigma2.^2).^(-1)) ;
    f2 = sum(sum(P1.*polamb2.*Phi2)) ;
    g = (f2-f1)./((sigma2-sigma1)*grad.') ;
    if g < eps1
        muk = muk/delta ;
        sigma2 = sigma1 ;
    elseif g > eps2
        muk = min(delta*muk,M) ;
    else
        criterion1 = norm(sigma2-sigma1)-sqrt(eps3)*(1+gamma) ;
        criterion2 = (f2-f1)-eps3*(1+f1) ;
    end
    n_iter = n_iter + 1 ;
end

%% Convert the RGD Kernel from polar to rectangular coordinates
psi = [psi psi(2:Q)+pi];
sigma = [sigma2 sigma2(2:Q)];
A = atan2(TAU,THETA);
negangle = find(A < 0);
A(negangle) = A(negangle)+2*pi;
R = TAU.^2+THETA.^2;
S = interp1(psi,sigma,A(:));
S = reshape(S,N,N);
G = exp(-R.*((2*S.^2).^(-1))) ;
G = G./max(max(G));
end