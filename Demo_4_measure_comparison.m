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
% This demo script produces the results that are depicted in Figs. 5 and 6
% of the paper. It compares the proposed average TFD performance with the
% NIR and Reinhold measures using the CKD of the example signal generated
% in Demo_2_tfd_comparison.m. In addition, it illustrates the TFD
% performance progression by computing the CKD at four increasing
% performance levels. Note that the CKD optimization is precomputed.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
N  = 256;
M  = 256;
fs = 1;
t  = 0:1/fs:(N-1)/fs;
f  = 0:fs/(2*N-1):fs/2;

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

%% PW-WVD (Lz)
Lz = pw_wvd(IF,IP,IA,W,fs);

%% WVD (Wz)
Wz = qtfd_wvd(z,N-1);

%% Compute Xz(t,f) and normlaize it by volume
Xz = Wz - Lz;
Xz = (sum(sum(abs(Xz))) < 1e-3).*Xz + (sum(sum(abs(Xz))) > 1e-3).*Xz./sum(sum(abs(Xz)));

%% Optimization (Uncomment the line below if you want to re-optimize the CKD)
% rng(1); rf = optimize_ckd_example(z,Lz,200,1);
rf = [2.99504082412799,0.395077450628600,0.424005135399697];

%% CKD (C,D) parameter range between r0 and rf
K  = 128;
p0 = 1;
re = linspace(p0,rf(3),K);
r0 = [rf(1), p0, p0];
m  = (rf(2)-r0(2))/(rf(3)-r0(3));
rd = m*re + rf(2) - m*rf(3);
r = [repmat(rf(1),K,1) rd' re'];

%% Calculate the proposed, NIR, and Reinhold measures
Pa   = zeros(1,K);
Pr   = zeros(1,K);
NIR  = zeros(1,K);
RNIR = zeros(1,K);
for k = 1:K
    G   = tf_kernel_ckd(N,r(k,1),r(k,2),r(k,3));
    TFa = filter_tfd(Xz,G);
    TFr = filter_tfd(Lz,G);
    TFD = filter_tfd(Wz,G);
    Pa(k)   = tfd_accuracy(TFa);
    Pr(k)   = tfd_resolution(Lz,TFr);
    NIR(k)  = tfd_measure(TFD');
    % The Reinhold measure does not work with TFDs that have low values;
    % hence, we excluded the first and last 10 samples in the computation
    RNIR(k) = adatnirperformance(TFD(10:N-10,:)');
end
P = 0.5.*(Pa + Pr);

%% Select the points of interest
ps = [40 80 K-9 K];
rs = r(ps,:);
TFD = cell(1,length(rs));
for i = 1:length(rs)
    G = tf_kernel(z,'ckd',rs(i,1),rs(i,2),rs(i,3));
    TFD{i} = filter_tfd(Wz,G);
end
disp('CKD parameters are:');
disp(r(ps,:));
disp('Proposed measure');
disp(round(P(ps),3));
disp('NIR measure');
disp(round(NIR(ps),3));
disp('Reinhold measure');
disp(round(RNIR(ps),3));

%% Plotting the proposed, NIR, and Reinhold measures
figure('Color',[1,1,1],'Position',[100 100 650 550]);
hAx(1) = gca;
plot(P,'-b','linewidth',3); hold on;
plot(NIR,'--k','linewidth',3); grid on;
plot(RNIR,'-.r','linewidth',3);
ylabel('Performance'); axis([1 K+1 0.5 0.9]);
hAx(1).OuterPosition = [0 0 1 0.93];
xlabel('$$D$$','Interpreter','LaTeX');
xticks(20:20:120); xticklabels(num2str(round(r(20:20:120,2),2)));
legend('$$\xi$$','NIR','Reinhold','Interpreter','LaTeX',...
    'location','southeast','fontsize',20,'FontName','Times',...
    'Orientation','horizontal');
set(gca,'fontweight','bold','fontsize',20,'FontName','Times');
hAx(1) = gca;
hAx(2) = axes('Position',hAx(1).Position,'XAxisLocation','top','Color','none');
hold(hAx(2),'on'); yticklabels('');
plot(ps,P(ps),'.b','Markersize',35);
plot(ps,NIR(ps),'.k','Markersize',35);
plot(ps,RNIR(ps),'.r','Markersize',35);
axis([1 K+1 0.5 0.9]); xlabel('$$E$$','Interpreter','LaTeX');
xticks(20:20:120); xticklabels(num2str(round(r(20:20:120,3),2)));
set(gca,'fontweight','bold','fontsize',20,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting the points of interest TFDs
C_max = 100;
for i = 1:length(TFD)
    figure('Color',[1,1,1],'Position',[100 100 650 550]);
    colormap(1-gray); imagesc(f,t,abs(TFD{i})); axis xy;
    xlabel('Frequency (Hz)'); ylabel('Time (s)');
    axis([f(1) f(end) t(1) t(end)]); caxis([0 C_max]);
    set(gca,'fontweight','bold','fontsize',24,'FontName','Times');
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
end

%% Saving
opt = input('Do you want to save all results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'measure_performance','-dpdf','-r400');
    print(2,'measure_TFD1','-dpdf','-r400');
    print(3,'measure_TFD2','-dpdf','-r400');
    print(4,'measure_TFD3','-dpdf','-r400');
    print(5,'measure_TFD4','-dpdf','-r400');
end