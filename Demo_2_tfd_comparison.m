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
% This demo script produces the results that are depicted in Fig. 3 of the
% paper. It generates an example 2-component non-stationary signal and
% computes the signal ideal TFD, PW-WVD, and the MPD using chirplet atoms.
% Note that the MPD is precomputed and saved in MPD_example.mat.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
N  = 256;
fs = 1;
t  = 0:1/fs:(N-1)/fs;
f  = 0:fs/(2*N-1):fs/2;
C_max = 128;

%% Generate component 1
fc = 0.3;
fd = 0.1;
fm = 0.015;
IF(1,:) = fc - fd.*sin(2*pi*fm*t);
IP(1,:) = fc*t + (fd./(2*pi*fm)).*cos(2*pi*fm*t);
W(1,:)  = ones(1,N);
IA(1,:) = ones(1,N);
z(1,:)  = IA(1,:).*exp(1j*2*pi*IP(1,:)).*W(1,:);

%% Generate component 2
ff = 0.025;
m  = 0.1/255;
IF(2,:) = ff + m*t;
IP(2,:) = ff*t + 0.5.*m.*t.^2;
W(2,:)  = ones(1,N);
IA(2,:) = exp(-0.0037*t);
z(2,:)  = IA(2,:).*exp(1j*2*pi*IP(2,:)).*W(2,:);
z       = sum(z);

%% Ideal TFD
Iz = ideal_tfd(IF,IA,W);

%% PW-WVD
Lz = pw_wvd(IF,IP,IA,W,fs);

%% MPD using chirplet atoms.
% Uncomment the 5 lines below if you want to re-compute the MPD.
% dict = chirplet_dictionary(N,fs);
% [atoms,c,idx] = mp_atoms(z,dict,1);
% Mz = mpd(atoms);
% save('Data\MP\MPD_example.mat','Mz');
% save('Data\MP\MP_example.mat','c','idx');
load('Data\MP\MPD_example.mat');

%% Ideal TFD TF-Tiling
ideal_t_tiling = repmat(t,N,1);

%% PW-WVD TF-Tiling
P = size(IF,1);
pw_wvd_tf_tiling = cell(1,P);
for p = 1:P
    IF_temp = IF(p,:);
    IFw = IF_temp - max(IF_temp);
    fd = 0:fs/(2*N-1):(abs(min(IFw))+fs/2);
    for i = 1:length(fd)
        pw_wvd_tf_tiling{p}(i,:) = IFw + fd(i);
    end
end

%% MPD TF-Tiling
load('Data\MP\MP_example.mat');
parameters = chirplet_parameters(N,fs);
mp_param = parameters(:,idx);
index = [1 4 8];
mpd_tf_tiling = cell(1,P);
for p = 1:length(index)
    IF_temp = mp_param(end,index(p)).*(t - mp_param(2,index(p)))./pi;
    IFw = IF_temp - max(IF_temp);
    fd = 0:fs/(2*N-1):(abs(min(IFw))+fs/2);
    for i = 1:length(fd)
        mpd_tf_tiling{p}(i,:) = IFw + fd(i);
    end
end

%% Plotting ideal TFD, PW-WVD, and MPD
x1 = 0.035; y1 = 40.5;  w1 = 0.022; h1 = 30;
x2 = 0.351;  y2 = 105.5; w2 = 0.057;  h2 = 22;

figure('Color',[1,1,1],'Position',[100 100 2200 550]);
f1 = subplot(1,3,1);
f1.Position = [0.125 0.1589 0.2134 0.7661];
colormap(1-gray); imagesc(f,t,abs(Iz)); axis xy; hold on;
plot(f,0.5+ideal_t_tiling(:,4:5:end),'-','linewidth',3,'Color',[0 0.4470 0.7410]);
rectangle('Position',[x1 y1 w1 h1],'edgecolor',[0 0.4470 0.7410],'linewidth',5);
axis([x1 x1+w1 y1 y1+h1]); caxis([0 C_max]);
set(gca,'XtickLabel','','YtickLabel','');
annotation('line',[0.338181818181818 0.425909090909091],...
    [0.923636363636363 0.367272727272727],'Color',[0 0.4470 0.7410],'LineWidth',5);
annotation('line',[0.338181818181818 0.425909090909091],[0.16 0.28],...
    'Color',[0 0.4470 0.7410],'LineWidth',5);
f2 = subplot(1,3,2);
colormap(1-gray); imagesc(f,t,abs(Iz)); axis xy; hold on;
rectangle('Position',[x1 y1 w1 h1],'edgecolor',[0 0.4470 0.7410],'linewidth',5);
rectangle('Position',[x2 y2 w2 h2],'edgecolor',[0.4660 0.6740 0.1880],'linewidth',5);
xlabel('Frequency (Hz)'); ylabel('Time (s)');
axis([f(1) f(end) t(1) t(end)]); caxis([0 C_max]);
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
f3 = subplot(1,3,3);
f3.Position = [0.7 0.1589 0.2134 0.7661];
colormap(1-gray); imagesc(f,t,abs(Iz)); axis xy; hold on;
plot(f,0.5+ideal_t_tiling,'-','linewidth',3,'Color',[0.4660 0.6740 0.1880]);
rectangle('Position',[x2 y2 w2 h2],'edgecolor',[0.4660 0.6740 0.1880],'linewidth',5);
axis([x2 x2+w2 y2 y2+h2]); caxis([0 C_max]);
set(gca,'XtickLabel','','YtickLabel','');
annotation('line',[0.585454545454545 0.7],...
    [0.543636363636364 0.925454545454545],'Color',[0.4660 0.6740 0.1880],'LineWidth',5);
annotation('line',[0.585 0.700454545454545],...
    [0.476363636363636 0.158181818181818],'Color',[0.4660 0.6740 0.1880],'LineWidth',5);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[100 100 2200 550]);
f1 = subplot(1,3,1);
f1.Position = [0.125 0.1589 0.2134 0.7661];
colormap(1-gray); imagesc(f,t,abs(Lz)); axis xy; hold on;
plot(0.5.*fs/(2*N-1)+pw_wvd_tf_tiling{2},t,'-','linewidth',3,'Color',[0 0.4470 0.7410]);
rectangle('Position',[x1 y1 w1 h1],'edgecolor',[0 0.4470 0.7410],'linewidth',5);
axis([x1 x1+w1 y1 y1+h1]); caxis([0 2*C_max]);
set(gca,'XtickLabel','','YtickLabel','');
annotation('line',[0.338181818181818 0.425909090909091],...
    [0.923636363636363 0.367272727272727],'Color',[0 0.4470 0.7410],'LineWidth',5);
annotation('line',[0.338181818181818 0.425909090909091],[0.16 0.28],...
    'Color',[0 0.4470 0.7410],'LineWidth',5);
f2 = subplot(1,3,2);
colormap(1-gray); imagesc(f,t,abs(Lz)); axis xy; hold on;
rectangle('Position',[x1 y1 w1 h1],'edgecolor',[0 0.4470 0.7410],'linewidth',5);
rectangle('Position',[x2 y2 w2 h2],'edgecolor',[0.4660 0.6740 0.1880],'linewidth',5);
xlabel('Frequency (Hz)'); ylabel('Time (s)');
axis([f(1) f(end) t(1) t(end)]); caxis([0 2*C_max]);
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
f3 = subplot(1,3,3);
f3.Position = [0.7 0.1589 0.2134 0.7661];
colormap(1-gray); imagesc(f,t,abs(Lz)); axis xy; hold on;
plot(0.5.*fs/(2*N-1)+pw_wvd_tf_tiling{1}(1:2:end,:),t,'-','linewidth',3,'Color',[0.4660 0.6740 0.1880]);
rectangle('Position',[x2 y2 w2 h2],'edgecolor',[0.4660 0.6740 0.1880],'linewidth',5);
axis([x2 x2+w2 y2 y2+h2]); caxis([0 2*C_max]);
set(gca,'XtickLabel','','YtickLabel','');
annotation('line',[0.585454545454545 0.7],...
    [0.543636363636364 0.925454545454545],'Color',[0.4660 0.6740 0.1880],'LineWidth',5);
annotation('line',[0.585 0.700454545454545],...
    [0.476363636363636 0.158181818181818],'Color',[0.4660 0.6740 0.1880],'LineWidth',5);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[100 100 2200 550]);
f1 = subplot(1,3,1);
f1.Position = [0.125 0.1589 0.2134 0.7661];
colormap(1-gray); imagesc(f,t,abs(Mz)); axis xy; hold on;
plot(0.5.*fs/(2*N-1)+mpd_tf_tiling{1},t,'-','linewidth',3,'Color',[0 0.4470 0.7410]);
rectangle('Position',[x1+0.002 y1 w1-0.002 h1],'edgecolor',[0 0.4470 0.7410],'linewidth',5);
axis([x1+0.002 x1+w1 y1 y1+h1]); caxis([0 C_max]);
set(gca,'XtickLabel','','YtickLabel','');
annotation('line',[0.338181818181818 0.425909090909091],...
    [0.923636363636363 0.367272727272727],'Color',[0 0.4470 0.7410],'LineWidth',5);
annotation('line',[0.338181818181818 0.425909090909091],[0.16 0.28],...
    'Color',[0 0.4470 0.7410],'LineWidth',5);
f2 = subplot(1,3,2);
colormap(1-gray); imagesc(f,t,abs(Mz)); axis xy; hold on;
rectangle('Position',[x1 y1 w1 h1],'edgecolor',[0 0.4470 0.7410],'linewidth',5);
rectangle('Position',[x2 y2 w2 h2],'edgecolor',[0.4660 0.6740 0.1880],'linewidth',5);
xlabel('Frequency (Hz)'); ylabel('Time (s)');
axis([f(1) f(end) t(1) t(end)]); caxis([0 C_max]);
set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
f3 = subplot(1,3,3);
f3.Position = [0.7 0.1589 0.2134 0.7661];
colormap(1-gray); imagesc(f,t,abs(Mz)); axis xy; hold on;
plot(0.5.*fs/(2*N-1)+mpd_tf_tiling{2}(1:5:end,1:117),t(1:117),'-','linewidth',3,'Color',[0.4660 0.6740 0.1880]);
plot(0.004+mpd_tf_tiling{3}(1:5:end,117:end),t(117:end),'-','linewidth',3,'Color',[0.4660 0.6740 0.1880]);
rectangle('Position',[x2 y2 w2 h2],'edgecolor',[0.4660 0.6740 0.1880],'linewidth',5);
axis([x2 x2+w2 y2 y2+h2]); caxis([0 0.3.*C_max]);
set(gca,'XtickLabel','','YtickLabel','');
annotation('line',[0.585454545454545 0.7],...
    [0.543636363636364 0.925454545454545],'Color',[0.4660 0.6740 0.1880],'LineWidth',5);
annotation('line',[0.585 0.700454545454545],...
    [0.476363636363636 0.158181818181818],'Color',[0.4660 0.6740 0.1880],'LineWidth',5);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save all results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'comparison_Iz','-dpdf','-r400');
    print(2,'comparison_Lz','-dpdf','-r400');
    print(3,'comparison_Mz','-dpdf','-r400');
end