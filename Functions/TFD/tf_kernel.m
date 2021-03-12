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
%                        TFDs Doppler-lag kernels
%
%  Syntax : G = tf_kernel(s,tfd,varargin)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% s        : The input signal (1 x N).
% tfd      : The TFD selection variable.
% varargin : The TFD kernel parameters.
%
% <OUTPUTs>
% G : The TFD Doppler-lag kernel (N x N).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function G = tf_kernel(s,tfd,varargin)
%% Generate TF kernel
N  = length(s);
switch lower(tfd)
    case 'wvd'
        G = ones(N,N);
    case 'pwvd'
        G = tf_kernel_pwvd(N,varargin{:});
    case 'spwvd'
        G = tf_kernel_spwvd(N,varargin{:});
    case 'ed'
        G = tf_kernel_ed(N,varargin{:});
    case 'bjd'
        G = tf_kernel_bjd(N,varargin{:});
    case 'bd'
        G = tf_kernel_bd(N,varargin{:});
    case 'mbd'
        G = tf_kernel_mbd(N,varargin{:});
    case 'embd'
        G = tf_kernel_embd(N,varargin{:});
    case 'ckd'
        G = tf_kernel_ckd(N,varargin{:});
    case 'rgd'
        G = tf_kernel_rgd(N,tf2af(qtfd_wvd(s,N-1)),varargin{:});
    case 'mdd'
        G = tf_kernel_mdd(N,tf2af(qtfd_wvd(s,N-1)),varargin{:});
end
end