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
% Last Modification: 21-May-2021
%
% Description:
% This demo script produces the results that are depicted in Fig. 8 of the
% paper. It produces the TFD computational complexity measured in terms of
% averaged processing time. It uses the 12 TFD computational complexity
% evaluations saved in comtime_pwvd.mat, comtime_spwvd.mat, comtime_ed.mat,
% comtime_bjd.mat, comtime_bd.mat, comtime_mbd.mat, comtime_embd.mat,
% comtime_ckd.mat, comtime_rgd.mat, comtime_mdd.mat, comtime_dgf.mat, and
% comtime_mpd.mat.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Main
tfd = {'pwvd','spwvd','ed','bjd','bd','mbd','embd','ckd','rgd','mdd','dgf','mpd'};
Tm = zeros(1,length(tfd));
Ts = zeros(1,length(tfd));
for i = 1:length(tfd)
    idx = tfd{i};
    load(['Data/Computation time/comtime_' idx '.mat']);
    Tm(i) = mean(T(:));
    Ts(i) = std(T(:));
end
[~,I] = sort(Tm,'descend');

%% Plotting
figure('Color',[1,1,1],'Position',[100 100 650 550]);
bar(Tm(I)); hold on;
errorbar(Tm(I),2.*Ts(I),'LineStyle','None','color','r',...
    'CapSize',12,'linewidth',1); grid on;
xlim([0.4 length(tfd)+0.5]); ylabel('Processing Time (s)')
xticks(1:length(tfd)); xticklabels(upper(tfd(I))); xtickangle(45);
set(gca,'fontweight','bold','fontsize',20,'FontName','Times','YScale','log');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save all results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'computation_time','-dpdf','-r400');
end