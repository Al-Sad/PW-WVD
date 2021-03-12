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
%                    Bayesian optimization variables
%
%  Syntax : [c,w] = prepare_opt_var(tfd,N)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% tfd : The TFD selection variable.
% N   : The signal total number of time samples.
%
% <OUTPUTs>
% c : Optimization parameters indices.
% w : Optimization parameters values.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [c,w] = prepare_opt_var(tfd,N)
L1  = 1:2:N-1;
win = {'rectwin','hanning','hamming','bartlett'};
L2  = 2:2:N-1;
switch lower(tfd)
    case 'pwvd'
        c(1) = optimizableVariable('L',[1 length(L1)],'Type','integer');
        c(2) = optimizableVariable('W',[1 length(win)],'Type','integer');
        w    = {L1, win};
    case 'spwvd'
        c(1) = optimizableVariable('L1',[1 length(L1)],'Type','integer');
        c(2) = optimizableVariable('W1',[1 length(win)],'Type','integer');
        c(3) = optimizableVariable('L2',[1 length(L1)],'Type','integer');
        c(4) = optimizableVariable('W2',[1 length(win)],'Type','integer');
        w    = {L1, win, L1, win};
    case 'spec'
        c(1) = optimizableVariable('L',[1 length(L1)],'Type','integer');
        c(2) = optimizableVariable('W',[1 length(win)],'Type','integer');
        w    = {L1, win};
    case 'ed'
        c(1) = optimizableVariable('alpha',[1e-6 1],'Type','real');
        w    = {};
    case 'bjd'
        c(1) = optimizableVariable('alpha',[0 1],'Type','real');
        w    = {};
    case 'bd'
        c(1) = optimizableVariable('alpha',[0 1],'Type','real');
        w    = {};
    case 'mbd'
        c(1) = optimizableVariable('alpha',[0 1],'Type','real');
        w    = {};
    case 'embd'
        c(1) = optimizableVariable('alpha',[0 1],'Type','real');
        c(2) = optimizableVariable('beta',[0 1],'Type','real');
        w    = {};
    case 'ckd'
        c(1) = optimizableVariable('c',[0 3],'Type','real');
        c(2) = optimizableVariable('d',[0.01 1],'Type','real');
        c(3) = optimizableVariable('e',[0.01 1],'Type','real');
        w    = {};
    case 'rgd'
        c = optimizableVariable('alpha',[1 20],'Type','real');
        w = {};
    case 'mdd'
        c(1) = optimizableVariable('c',[0 3],'Type','real');
        c(2) = optimizableVariable('thr',[0 1],'Type','real');
        w    = {};
    case 'aok'
        c = optimizableVariable('vol',[1 20],'Type','real');
        w = {};
    case 'dgf'
        c(1) = optimizableVariable('alpha1',[0 15],'Type','real');
        c(2) = optimizableVariable('alpha2',[0 15],'Type','real');
        c(3) = optimizableVariable('NNL',[1 length(L2)],'Type','integer');
        w    = {L2};
end
end