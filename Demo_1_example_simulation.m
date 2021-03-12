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
% This demo script produces the results that are depicted in Fig. 2 of the
% paper. It generates an example 2-component non-stationary signal and
% computes its PW-WVD, WVD and the WVD cross-terms. In addition, it shows
% the signal IA, IF, and time-support functions.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
N  = 256;
M  = 256;
fs = 1;
t  = 0:1/fs:(N-1)/fs;
f  = 0:fs/(2*N-1):fs/2;
C_max = 50;

%% Generate component 1
pt  = [20 80 130 180 230]';
pf  = [0.4 0.2 0.35 0.15 0.3]';
w   = (t >= min(pt) & t <= max(pt));
R   = length(pt) - 1;
tif = zeros(N,R+1);
tif(:,end) = ones(N,1);
for r = 1:R
    tif(:,r) = t.^(R-r+1);
end
tip = [(t.^(R+1))' tif(:,1:R)]./((R+1):-1:1);
RR = repmat(R:-1:0,R+1,1);
A  = pt.^RR;
C  = A\pf;
IF(1,:) = (tif*C)';
IP(1,:) = (tip*C)';
IA(1,:) = exp(-0.008*t);
W(1,:)  = w;
z(1,:)  = IA(1,:).*exp(1j*2*pi*IP(1,:)).*W(1,:);

%% Generate component 2
pt  = [100 220]';
pf  = [0.05 0.4]';
w   = (t >= min(pt) & t <= max(pt));
R   = length(pt) - 1;
tif = zeros(N,R+1);
tif(:,end) = ones(N,1);
for r = 1:R
    tif(:,r) = t.^(R-r+1);
end
tip = [(t.^(R+1))' tif(:,1:R)]./((R+1):-1:1);
RR  = repmat(R:-1:0,R+1,1);
A   = pt.^RR;
C   = A\pf;
IF(2,:) = (tif*C)';
IP(2,:) = (tip*C)';
IA(2,:) = 0.4 + 0.3.*cos(0.01*pi*t);
W(2,:)  = w;
z(2,:)  = IA(2,:).*exp(1j*2*pi*IP(2,:)).*W(2,:);

%% PW-WVD (Lz)
Lz = pw_wvd(IF,IP,IA,W,fs);

%% WVD (Wz)
Wz = qtfd_wvd(sum(z),N-1);

%% WVD cross-terms (Xz)
Xz = Lz - Wz;

%% Plotting IA laws
figure('Color',[1,1,1],'Position',[100 100 650 550]);
plot(t,IA(1,:),'b','linewidth',3); hold on;
plot(t,IA(2,:),'r','linewidth',3); hold on;
legend('$$a_1(t)$$','$$a_2(t)$$','Interpreter','LaTeX','location',...
    'northeast','fontsize',30,'FontName','Times','orientation','horizontal',...
    'Position',[0.18420512453715 0.827878790479718 0.411382203954988 0.0945454519445246]);
grid on; axis([t(1) t(end) 0 1.05]);
xlabel('Time (s)'); ylabel('Amplitude');
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting IF laws
figure('Color',[1,1,1],'Position',[100 100 650 550]);
plot(IF(1,:),t,'b','linewidth',3); hold on;
plot(IF(2,:),t,'r','linewidth',3); hold on;
legend('$$f_1(t)$$','$$f_2(t)$$','Interpreter','LaTeX','location',...
    'northeast','fontsize',30,'FontName','Times','orientation','horizontal',...
    'Position',[0.18420512453715 0.827878790479718 0.411382203954988 0.0945454519445246]);
grid on; axis([f(1) f(end) t(1) t(end)]);
xlabel('Frequency (Hz)'); ylabel('Time (s)');
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting Windows
figure('Color',[1,1,1],'Position',[100 100 650 550]);
plot(t,W(1,:),'b','linewidth',3); hold on;
plot(t,W(2,:),'r','linewidth',3); hold on;
legend('$$\Pi_1(t)$$','$$\Pi_2(t)$$','Interpreter','LaTeX','location',...
    'northeast','fontsize',30,'FontName','Times','orientation','horizontal',...
    'Position',[0.18420512453715 0.827878790479718 0.411382203954988 0.0945454519445246]);
grid on; axis([t(1) t(end) -0.05 1.19]);
xlabel('Time (s)'); ylabel('Amplitude');
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting Lz
figure('Color',[1,1,1],'Position',[100 100 650 550]);
colormap(1-gray); imagesc(f,t,abs(Lz)); axis xy;
xlabel('Frequency (Hz)'); ylabel('Time (s)');
axis([f(1) f(end) t(1) t(end)]); caxis([0 C_max]);
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting Wz
figure('Color',[1,1,1],'Position',[100 100 650 550]);
colormap(1-gray); imagesc(f,t,abs(Wz)); axis xy;
xlabel('Frequency (Hz)'); ylabel('Time (s)');
axis([f(1) f(end) t(1) t(end)]); caxis([0 C_max]);
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting Xz
figure('Color',[1,1,1],'Position',[100 100 650 550]);
colormap(1-gray); imagesc(f,t,abs(Xz)); axis xy;
xlabel('Frequency (Hz)'); ylabel('Time (s)');
axis([f(1) f(end) t(1) t(end)]); caxis([0 C_max]);
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save all results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'IA','-dpdf','-r400');
    print(2,'IF','-dpdf','-r400');
    print(3,'IW','-dpdf','-r400');
    print(4,'Lz','-dpdf','-r400');
    print(5,'Wz','-dpdf','-r400');
    print(6,'Xz','-dpdf','-r400');
end