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
% Description:
% This demo script produces the results that are depicted in Fig. A.1 of
% the paper. It generates an example mono-component non-stationary signal
% with non-linear IF and computes the WVD auto-terms and inner-terms using
% the formulation in Appendix A.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
fs  = 1;
N   = 256;
M   = N;
t   = 0:1/fs:(N-1)/fs;
tau = linspace(0,2*N,N+1);
f   = 0:fs/(2*N-1):fs/2;
C_max = 64;

%% Generate the signal polynomial IF
tt  = [t(1) t(end/2) t(end)];
ff  = [0.2 0.4 0.2];
R   = length(tt)-1;
P   = @(x) [x^2 x 1];
A   = [P(tt(1)) ; P(tt(2)) ; P(tt(3))];
eta = flipud(A\(ff'));
IF  = 0;
for r = 1:R+1
    IF  = IF + eta(r).*t.^(r-1);
end

%% Compute the recursive residual term zeta(t,\tau)
zeta_R = 0;
for r = 2:R+1
    for i = 2:r
        zeta_R = zeta_R + (eta(r).*(t.^(r-i))./r)'.*nchoosek(r,i)*((tau./(2)).^i).*(2.*rem(i,2));
    end
end

%% Compute K_pi(t,\tau)
K_pi = zeros(N,N);
for i = 0:(N/2)
    K_pi(:,i+1) = [zeros(i,1) ; ones(N-2*i,1) ;zeros(i,1)];
end
K_pi = [K_pi(:,1:N/2) zeros(N,1) fliplr(K_pi(:,2:N/2))];

%% Compute the WVD auto-terms W_a(t,f) * \delta(t,f-f(t)) * W_pi(t,f)
K_f   = exp(1j*2*pi*IF'*tau);
Phi_f = K_pi.*[K_f(:,1:N/2) zeros(N,1) fliplr(real(K_f(:,2:N/2))-1j*imag(K_f(:,2:N/2)))];
Phi_f = [Phi_f zeros(N,N); zeros(N,2*N)];
TF_f  = real(fft(Phi_f,N,2));
WVD_f = TF_f(1:N,:);

%% Compute the WVD inner-terms W_a(t,f) * W_zeta(t,f) * W_pi(t,f)
K_zeta = exp(1j*2*pi*zeta_R);
Phi_g  = K_pi.*[K_zeta(:,1:N/2) zeros(N,1) fliplr(real(K_zeta(:,2:N/2))-1j*imag(K_zeta(:,2:N/2)))];
Phi_g  = [Phi_g zeros(N,N); zeros(N,2*N)];
TF_g   = real(fftshift(fft(Phi_g,N,2),2));
WVD_g  = TF_g(1:N,:);

%% Compute the complete WVD 
K_c  = K_f.*K_zeta;
Phi_z = K_pi.*[K_c(:,1:N/2) zeros(N,1) fliplr(real(K_c(:,2:N/2))-1j*imag(K_c(:,2:N/2)))];
Phi_z = [Phi_z zeros(N,N); zeros(N,2*N)];
TF_z  = real(fft(Phi_z,N,2));
WVD_z = TF_z(1:N,:);

%% Plotting IF law
figure('Color',[1,1,1],'Position',[100 100 650 550]);
plot(IF,t,'b','linewidth',3); grid on;
axis([f(1) f(end) t(1) t(end)]);
xlabel('Frequency (Hz)'); ylabel('Time (s)');
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting the WVD auto-terms
figure('Color',[1,1,1],'Position',[100 100 650 550]);
colormap(1-gray); imagesc(f,t,abs(WVD_f)); axis xy;
xlabel('Frequency (Hz)'); ylabel('Time (s)');
axis([f(1) f(end) t(1) t(end)]); caxis([0 C_max]);
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting the WVD inner-terms
figure('Color',[1,1,1],'Position',[100 100 650 550]);
colormap(1-gray); imagesc(-fs/2:fs/(N-1):fs/2,t,abs(WVD_g)); axis xy;
xlabel('Frequency (Hz)'); ylabel('Time (s)');
axis([-fs/2 fs/2 t(1) t(end)]); caxis([0 C_max]);
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting the complete WVD
figure('Color',[1,1,1],'Position',[100 100 650 550]);
colormap(1-gray); imagesc(f,t,abs(WVD_z)); axis xy;
xlabel('Frequency (Hz)'); ylabel('Time (s)');
axis([f(1) f(end) t(1) t(end)]); caxis([0 C_max]);
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save all results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'appendix_IF','-dpdf','-r400');
    print(2,'appendix_auto_terms','-dpdf','-r400');
    print(3,'appendix_inner_terms','-dpdf','-r400');
    print(4,'appendix_Wz','-dpdf','-r400');
end