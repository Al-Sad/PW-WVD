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
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Bayesian optimization for the example CKD
%
%  Syntax : x = optimize_ckd_example(z,Lz,N_iter,use_parallel)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% z            : The input signal (1 x N).
% Lz           : The signal PW-WVD (N x N).
% N_iter       : Total number of iterations.
% use_parallel : 1 or 0 to perform parallel computations.
%
% <OUTPUTs>
% x : The CKD optimal parameters.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function x = optimize_ckd_example(z,Lz,N_iter,use_parallel)
%% WVD (Wz)
Wz = qtfd_wvd(z,length(z)-1);

%% Main
opt_param(1) = optimizableVariable('c',[0 3],'Type','real');
opt_param(2) = optimizableVariable('d',[0.01 1],'Type','real');
opt_param(3) = optimizableVariable('e',[0.01 1],'Type','real');
Initial_x    = table(1,1,1);
Fun = @(opt_param)objective_example_fun(z,Wz,Lz,opt_param);
bayes_out = bayesopt(Fun,opt_param,'Verbose',1,'UseParallel',use_parallel,...
    'PlotFcn',[],'MaxObjectiveEvaluations',N_iter,...
    'AcquisitionFunctionName','expected-improvement-plus',...
    'IsObjectiveDeterministic',0,'InitialX',Initial_x);
x = table2array(bayes_out.XAtMinObjective);

%% Cost function
function cost = objective_example_fun(z,Wz,Lz,opt_param)
G  = tf_kernel(z,'ckd',opt_param.c,opt_param.d,opt_param.e);
TF = filter_tfd(Wz,G);
cost = sum(sum((Lz - TF).^2));
end
end