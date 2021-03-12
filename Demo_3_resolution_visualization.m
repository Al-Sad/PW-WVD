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
% This demo script produces the illustration in Fig. 4 of the paper. It
% produces a one-dimensional visual interpretation for the proposed TFD
% resolution measure.

%% Initialization
clear; close all; clc;

%% Parameters
N  = 2^12;
M  = 2^12;
fs = 1;
t = linspace(-100,100,N);

%% Main
T = 10;
x = sinc(t/T).*(t>=-2*T & t<=2*T);
x = x./sum(abs(x));
y = smoothdata(x,'gaussian','SmoothingFactor',0.7);
z = sum(min(abs(x),abs(y)));
I = abs(x) <= abs(y);
zz = zeros(size(x));
zz(I)  = x(I);
zz(~I) = y(~I);
disp(round(z,3));

%% Plotting
figure('Color',[1,1,1],'Position',[100 100 1300 450]);
f1 = patch([t fliplr(t)], [zeros(1,N) min(abs(x),abs(y))],...
    [0.9,0.9,0.9],'linestyle','none'); hold on;
f2 = plot(t,abs(x),'linewidth',3); hold on;
f3 = plot(t,abs(y),'linewidth',3);
xticks({}); yticks({}); xticklabels({}); yticklabels({});
axis([-30 30 -1e-5 1.8e-4+max([x,y])]); box on;
legend([f2 f3 f1],'$$\Big|\frac{L_z(f)}{\int |L_z(f)| df}\Big|$$','$$|\ell_z(f)|$$','Intersection area',...
    'Interpreter','LaTeX','location','northeast',...
    'fontsize',28,'FontName','Times','Orientation','vertical');
xlabel('Frequency (Hz)');
set(gca,'fontweight','bold','fontsize',20,'FontName','Times');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save all results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'resolution_visual','-dpdf','-r400');
end