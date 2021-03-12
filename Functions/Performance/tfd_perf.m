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
%                    Proposed TFD performance measure
%
%  Syntax : [Pa, Pr] = tfd_perf(z,Lz,tfd,varargin)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% z        : The input signal (1 x N).
% Lz       : The PW-WVD (N x N).
% tfd      : The TFD selection variable.
% varargin : The TFD kernel parameters.
%
% <OUTPUTs>
% Pa : The TFD accuracy between 0 and 1.
% Pr : The TFD resolution between 0 and 1.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [Pa, Pr] = tfd_perf(z,Lz,tfd,varargin)
%% Compute Wz(t,f)
Wz = qtfd_wvd(z,length(z)-1);

%% Compute Xz(t,f) and normlaize it by volume
Xz = Wz - Lz;
Xz = (sum(sum(abs(Xz))) < 1e-3).*Xz + (sum(sum(abs(Xz))) > 1e-3).*Xz./sum(sum(abs(Xz)));

%% Compute optimal kernel
switch lower(tfd)
    case 'dgf'
        TFa = tf_kernel_dgf(Xz,Wz,varargin{:});
        TFr = tf_kernel_dgf(Lz,Wz,varargin{:});
        Pa  = tfd_accuracy(TFa);
        Pr  = tfd_resolution(Lz,TFr);
    case 'mpd'
        Pa  = 1;
        TFr = varargin{1};
        Pr  = tfd_resolution(Lz,TFr);
    otherwise
        G   = tf_kernel(z,tfd,varargin{:});
        TFa = filter_tfd(Xz,G);
        TFr = filter_tfd(Lz,G);
        Pa  = tfd_accuracy(TFa);
        Pr  = tfd_resolution(Lz,TFr);
end
end