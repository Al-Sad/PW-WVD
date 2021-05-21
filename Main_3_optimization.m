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
% This main script optimizes each quadratic TFD using the Bayesian
% optimization algorithm. The optimization is executed with random initial
% kernel parameters, for 200 iterations, and by using the expected
% improvement plus acquisition function. The optimal TFD parameters are
% then saved in opt_pwvd.mat, opt_spwvd.mat, opt_ed.mat, opt_bjd.mat,
% opt_bd.mat, opt_mbd.mat, opt_embd.mat, opt_ckd.mat, opt_rgd.mat,
% opt_mdd.mat, and opt_dgf.mat.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
N_iter  = 200;
Use_par = 1;
disp_re = 0;

%% Loading database
load(['Data' filesep 'Database' filesep 'pw_wvd_database.mat'],'z','N','K','IF','IP','IA','Wp','fs');

%% Compute the PW-WVD and WVD of each signal z
Lz = cell(1,K);
Wz = cell(1,K);
for k = 1:K
    Lz{k} = pw_wvd(IF{k},IP{k},IA{k},Wp{k},fs);
    Wz{k} = qtfd_wvd(z{k},N-1);
end

%% Main
tfd = {'pwvd','ed','bjd','bd','mbd','embd','ckd','rgd','mdd','dgf','aok'};
for i = 1:length(tfd)
    idx = tfd{i};
    fprintf('%s \n',idx);
    [opt_var,w] = prepare_opt_var(idx,N);
    param = cell(K,length(opt_var));
    for k = 1:K
        fprintf('%0.2f %%\n',100*k/K);
        rng(1);
        Fun = @(opt_var)objective_fun(z{k},Wz{k},Lz{k},idx,w,opt_var);
        bayes_out = bayesopt(Fun,opt_var,'Verbose',disp_re,'UseParallel',Use_par,...
            'PlotFcn',[],'MaxObjectiveEvaluations',N_iter,...
            'AcquisitionFunctionName','expected-improvement-plus',...
            'IsObjectiveDeterministic',0);
        param(k,:) = tfd_param(idx,bayes_out.XAtMinObjective,w);
    end
    save(['Data' filesep 'Optimization' filesep 'opt_' idx '.mat'],'param');
end

%% Optimize SPWVD using PWVD optimal parameters as initial guesses
idx = 'spwvd';
load(['Data' filesep 'Optimization' filesep 'opt_pwvd.mat']);
param_pwvd = param; clear param;
[opt_var,w] = prepare_opt_var(idx,N);
I = zeros(K,2);
for k = 1:K
    I(k,1) = find((w{1} == param_pwvd{k,1}));
    I(k,2) = find(strcmp(param_pwvd{k,2},{'rectwin','hanning','hamming','bartlett'}));
end
Initial_x = table(I(:,1),I(:,2),ones(K,1),ones(K,1));
param = cell(K,length(opt_var));
for k = 1:K
    fprintf('%0.2f %%\n',100*k/K);
    rng(1);
    Fun = @(opt_var)objective_fun(z{k},Wz{k},Lz{k},idx,w,opt_var);
    bayes_out = bayesopt(Fun,opt_var,'Verbose',disp_re,'UseParallel',Use_par,...
        'PlotFcn',[],'MaxObjectiveEvaluations',N_iter,...
        'AcquisitionFunctionName','expected-improvement-plus',...
        'IsObjectiveDeterministic',0,'InitialX',Initial_x(k,:));
    param(k,:) = tfd_param(idx,bayes_out.XAtMinObjective,w);
end
save(['Data' filesep 'Optimization' filesep 'opt_' idx '.mat'],'param');