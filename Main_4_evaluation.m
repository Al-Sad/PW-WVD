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
% This main script evaluates each optimized TFD using the proposed accuracy
% and resolution measures. The computed evaluations are saved in
% perf_pwvd.mat, perf_spwvd.mat, perf_ed.mat, perf_bjd.mat, perf_bd.mat,
% perf_mbd.mat, perf_embd.mat, perf_ckd.mat, perf_rgd.mat, perf_mdd.mat,
% perf_dgf.mat, and perf_mpd.mat.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Loading database
load('Data\Database\pw_wvd_database.mat','z','N','K','IF','IP','IA','Wp','fs');

%% Compute the PW-WVD of each signal z
Lz = cell(1,K);
for k = 1:K
    disp(100*k/K)
    Lz{k} = pw_wvd(IF{k},IP{k},IA{k},Wp{k},fs);
end

%% Main
tfd = {'pwvd','spwvd','ed','bjd','bd','mbd','embd','ckd','rgd','mdd','dgf','mpd'};
for i = 1:length(tfd)
    idx = tfd{i};
    fprintf('%s \n',idx);
    Pa = zeros(1,K);
    Pr = zeros(1,K);
    if(~strcmp(idx,'mpd'))
        load(['Data\Optimization\opt_' idx '.mat']);
        for k = 1:K
            disp(k/K*100)
            [Pa(k), Pr(k)] = tfd_perf(z{k},Lz{k},idx,param{k,:});
        end
    else
        load('Data\MP\MP_atoms.mat');
        for k = 1:K
            disp(k/K*100)
            TFD = mpd(Atoms{k});
            [Pa(k), Pr(k)] = tfd_perf(z{k},Lz{k},idx,TFD);
        end
    end
    save(['Data\Evaluation\perf_' idx '.mat'],'Pa','Pr');
end
