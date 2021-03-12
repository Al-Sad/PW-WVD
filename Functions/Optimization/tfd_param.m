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
%           Getting TFDs parameters from the Bayesian Optimizer
%
%  Syntax : param = tfd_param(tfd,opt_var,w)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% tfd     : The TFD selection variable.
% opt_var : Optimization parameters indices.
% w       : Optimization parameters values.
%
% <OUTPUTs>
% param : The TFD parameters.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function param = tfd_param(tfd,opt_var,w)
switch lower(tfd)
    case 'pwvd'
        param{1} = w{1}(opt_var.L);
        param{2} = w{2}{opt_var.W};
    case 'spwvd'
        param{1} = w{1}(opt_var.L1);
        param{2} = w{2}{opt_var.W1};
        param{3} = w{3}(opt_var.L2);
        param{4} = w{4}{opt_var.W2};
    case 'ed'
        param{1} = opt_var.alpha;
    case 'bjd'
        param{1} = opt_var.alpha;
    case 'bd'
        param{1} = opt_var.alpha;
    case 'mbd'
        param{1} = opt_var.alpha;
    case 'embd'
        param{1} = opt_var.alpha;
        param{2} = opt_var.beta;
    case 'ckd'
        param{1} = opt_var.c;
        param{2} = opt_var.d;
        param{3} = opt_var.e;
    case 'rgd'
        param{1} = opt_var.alpha;
    case 'mdd'
        param{1} = opt_var.c;
        param{2} = opt_var.thr;
    case 'aok'
        param{1} = opt_var.vol;
    case 'dgf'
        param{1} = opt_var.alpha1;
        param{2} = opt_var.alpha2;
        param{3} = w{1}(opt_var.NNL);
end
end