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
% supplementary material of Boualem Boashash, Time-frequency signal
% analysis and processing toolbox, Online (2016). URL: http://booksite.
% elsevier.com/9780123984999/toolbox.php.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               The directional Gaussian filter distribution
%
%  Syntax : TF = qtfd_dgf(tf,alpha1,alpha2,N2)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% tf     : The input WVD (N x N).
% alpha1 : The smoothing parameter along the major ridge axis.
% alpha2 : The smoothing parameter along the minor ridge axis.
% N2     : An integer between 1 and N that defines the number of filtering
%          directions.
%
% <OUTPUTs>
% TF : The directional Gaussian filter distribution of the input
%      signal (N x N).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function TF = qtfd_dgf(tf,alpha1,alpha2,N2)
%% Parameters
N = size(tf,1);
N_theta = 60;

%% Pre-filtering using the MBD
alpha = (N2 <= 200)*0.6 + (N2 > 200)*0.25;
t  = -floor(N2/2):floor(N2/2);
h = 1./((cosh(t).^2).^alpha);
h = h/sum(h);
B_mbd = h'*h;
tf = filter2(B_mbd,tf);

%% Generate a set of directional filters
[X,Y] = meshgrid(-1:2/N2:1,-1:2/N2:1);
B = cell(N_theta,1);
for i = 1:N_theta
    angle = pi*(i-1)/N_theta;
    Xr = X*cos(angle)-Y*sin(angle);
    Yr = X*sin(angle)+Y*cos(angle);
    A  = exp((-1/2)*(((alpha1*Xr).^2)+(alpha2*Yr).^2));
    A  = A.*(1-alpha2*alpha2*Yr.^2);
    A  = A/sum(sum(abs(A)));
    B{i} = A;
end

%% Apply directional filters
M  = 2*N;
Q  = round((N2+1)/2);
F  = fft2(double(tf), M, M);
G  = fft2(double((abs(tf))), M, M);
I1 = zeros([M-length(A) M-length(A) N_theta]);
I2 = zeros([M-length(A) M-length(A) N_theta]);
for i = 1:N_theta 
    H = fft2(B{i}, M, M);
    F_fH = H.*F;
    G_fH = H.*G;
    tf_r = real(ifft2(F_fH));
    Wr   = (abs(ifft2(G_fH))).^2;
    I1(:,:,i)= Wr(Q:end-Q, Q:end-Q);
    I2(:,:,i) = tf_r(Q:end-Q, Q:end-Q);
end

%% Find the filter that gives the highest ampltiude
[~,a] = max(I1,[],3);

%% Computing the Directional Gaussian Filter Distribution
TF = tf;
for m = 1:N
    for n = 1:N
        TF(m,n) = I2(m,n,a(m,n));
    end
end
TF = filter2(B_mbd,TF);
end