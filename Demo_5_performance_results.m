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
% This demo script produces the results that are depicted in Fig. 7 and
% Table I of the paper. It yields the TFD performance evaluation results
% using the proposed measures. It uses the database signal parameters saved
% in pw_wvd_database.mat and the 12 TFD evaluations saved in perf_pwvd.mat,
% perf_spwvd.mat, perf_ed.mat, perf_bjd.mat, perf_bd.mat, perf_mbd.mat,
% perf_embd.mat, perf_ckd.mat, perf_rgd.mat, perf_mdd.mat, perf_dgf.mat,
% and perf_mpd.mat.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Loading database information
load('Data/Database/pw_wvd_database.mat','S','PR','K');

%% Main
tfd = {'pwvd','spwvd','ed','bjd','bd','mbd','embd','ckd','rgd','mdd','dgf','mpd'};
Perf_a = zeros(K,length(tfd));
Perf_r = zeros(K,length(tfd));
Perf_m = zeros(K,length(tfd));
for i = 1:length(tfd)
    idx = tfd{i};
    load(['Data/Evaluation/perf_' idx '.mat']);
    Perf_a(:,i) = Pa;
    Perf_r(:,i) = Pr;
    Perf_m(:,i) = (Pa + Pr)./2;
end
accuracy_md    = mean(Perf_a);
resolution_md  = mean(Perf_r);
performance_md = mean(Perf_m);
accuracy_sd    = std(Perf_a);
resolution_sd  = std(Perf_r);
performance_sd = std(Perf_m);
[~,I] = sort(performance_md,'ascend');

%% TFD selection results
mC = zeros(length(tfd),max(S));
for i = 1:length(tfd)
    for j = 1:max(S)
        mC(i,j) = mean(Perf_m(S == j,i));
    end
end
F = [];
for p = 1:size(PR,2)
    F = [F ; zeros(int8(nchoosek(p+max(max(PR))-1,p)),size(PR,2)-p) combs_rep(max(max(PR)),p)];
end
L = cell(1,size(F,1));
for i = 1:size(F,1)
    tt = [];
    Zn = sum(F(i,:)==0);
    for j = (Zn+1):size(PR,2)-1
        tt = [tt num2str(F(i,j)) ','];
    end
    tt = [tt num2str(F(i,size(PR,2)))];
    L{i} = ['(' tt ')'];
end

%% Find the set of best performing TFDs
thr = 0.01;
y = round(mC(I,:)',3);
c = repmat(max(y,[],2),1,length(tfd))-thr;
max_tfd = y >= c;

%% Latex results
header = '$\bm{P}$ & $\bm{R}$';
for i = 1:length(tfd)
    header = [header ' & \textbf{' upper(tfd{I(i)}) '}'];
end
header = [header ' \\'];
out = mC(I,:)';
OS = length(F);
output = cell(size(F,1),1);
% disp(header);
for j = 1:size(F,1)
    temp = ['& \textbf{' char(L{j}) '}'];
    for i = 1:length(tfd)
        if(max_tfd(j,i))
            temp = [temp ' & \cellcolor{Gray}\textbf{' num2str(round(out(j,i),3)) '}'];
        else
            temp = [temp ' & ' num2str(round(out(j,i),3))];
        end
    end
    output{j} = [temp ' \\ \hhline{|~|*{' num2str(length(tfd)+1) '}{-}}'];
%     disp(output{j})
end

%% Table I results
fprintf('\t\t\t\t %s',upper(tfd{I(1)}));
fprintf('\t  %s',upper(tfd{I(2)}));
fprintf('\t %s',upper(tfd{I(3)}));
fprintf('\t%s',upper(tfd{I(4)}));
fprintf('\t %s',upper(tfd{I(5)}));
fprintf('\t %s',upper(tfd{I(6)}));
fprintf('\t\t %s',upper(tfd{I(7)}));
fprintf('\t %s',upper(tfd{I(8)}));
fprintf('\t %s',upper(tfd{I(9)}));
fprintf('\t %s',upper(tfd{I(10)}));
fprintf('\t%s',upper(tfd{I(11)}));
fprintf('\t %s',upper(tfd{I(12)}));
fprintf('\n');
for j = 1:size(F,1)
    switch length(L{j})
        case 3, p = 1;
        case 5, p = 2;
        case 7, p = 3;
        case 9, p = 4;
    end
    ss = repelem(' ',1,12-length(L{j}));
    fprintf(['%d\t%s' ss],p,L{j});
    for i = 1:length(tfd)
        if(max_tfd(j,i))
            fprintf(2,'%0.3f\t',round(out(j,i),3));
        else
            fprintf('%0.3f\t',round(out(j,i),3));
        end
    end
    fprintf('\n');
end

%% Plotting
figure('Color',[1,1,1],'Position',[100 100 750 570]);
errorbar(accuracy_md(I), accuracy_sd(I),'LineStyle','None','color','b',...
    'CapSize',12,'linewidth',1); hold on;
p1 = plot(accuracy_md(I),'b-s','linewidth',3); hold on;
errorbar(resolution_md(I), resolution_sd(I),'LineStyle','None','color','r',...
    'CapSize',12,'linewidth',1); hold on;
p2 = plot(resolution_md(I),'r-o','linewidth',3); hold on;
p3 = plot(performance_md(I),'k-d','linewidth',3); grid on;
yy = accuracy_md(I);
text(1-0.2,yy(1)-0.03,num2str(round(yy(1),3)),'Interpreter','LaTeX',...
    'fontsize',12,'FontName','Times','color','b');
for i = [2,5,8,9]
    text(i+0.05,yy(i)+0.03,num2str(round(yy(i),3)),'Interpreter','LaTeX',...
        'fontsize',12,'FontName','Times','color','b');
end
for i = [3,4,6,7,10:12]
    text(i+0.05,yy(i)-0.04,num2str(round(yy(i),3)),'Interpreter','LaTeX',...
        'fontsize',12,'FontName','Times','color','b');
end
yy = resolution_md(I);
text(1+0.2,yy(1)+0.01,num2str(round(yy(1),3)),'Interpreter','LaTeX',...
    'fontsize',12,'FontName','Times','color','r');
for i = [4,7,11,12]
    text(i+0.05,yy(i)+0.03,num2str(round(yy(i),3)),'Interpreter','LaTeX',...
        'fontsize',12,'FontName','Times','color','r');
end
for i = [2,3,5,6,8:10]
    text(i+0.05,yy(i)-0.03,num2str(round(yy(i),3)),'Interpreter','LaTeX',...
        'fontsize',12,'FontName','Times','color','r');
end
yy = performance_md(I);
for i = 2:2:length(tfd)
    text(i,yy(i)+0.025,num2str(round(yy(i),3)),'Interpreter','LaTeX',...
        'fontsize',12,'FontName','Times','color','k');
end
for i = 1:2:length(tfd)
    text(i,yy(i)-0.03,num2str(round(yy(i),3)),'Interpreter','LaTeX',...
        'fontsize',12,'FontName','Times','color','k');
end
legend([p1 p2 p3],'$$\xi_a$$','$$\xi_r$$','$$\xi$$','Interpreter','LaTeX','location','southeast',...
    'fontsize',24,'FontName','Times','Orientation','horizontal');
axis([0.75 length(tfd)+1 0.175 1.01]); yticks(0:0.05:1); yticklabels(0:0.05:1);
xticks(1:length(tfd)); xticklabels(upper(tfd(I))); xtickangle(45);
set(gca,'fontweight','bold','fontsize',20,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save all results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'evaluation','-dpdf','-r400');
end