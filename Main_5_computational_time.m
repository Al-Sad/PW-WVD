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
% This main script evaluates the computational complexity of each optimized
% TFD in terms of processing time. The TFD computational complexity is
% estimated by Monte-Carlo simulations where the processing time of each
% optimized TFD is judged by generating the TFR of the 1000 test signals
% and repeating the process for 10 times. The TFD computational complexity
% evaluations are then saved in comtime_pwvd.mat, comtime_spwvd.mat,
% comtime_ed.mat, comtime_bjd.mat, comtime_bd.mat, comtime_mbd.mat,
% comtime_embd.mat, comtime_ckd.mat, comtime_rgd.mat, comtime_mdd.mat,
% comtime_dgf.mat, and comtime_mpd.mat.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Loading database
load('Data/Database/pw_wvd_database.mat','z','N','fs','K');

%% Main
L = 10; % Number of repetitions
tfd = {'pwvd','spwvd','ed','bjd','bd','mbd','embd','ckd','rgd','mdd','dgf','mpd'};
for i = 1:length(tfd)
    idx = tfd{i};
    fprintf('%s \n',idx);
    if(~strcmp(idx,'mpd'))
        load(['Data/Optimization/opt_' idx '.mat']);
    else
        dict = chirplet_dictionary(N,fs);
    end
    switch idx
        case 'pwvd'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    G  = tf_kernel_pwvd(N,param{k,:});
                    TF = filter_tfd(qtfd_wvd(z{k},N-1),G);
                    T(l,k) = toc;
                end
            end
        case 'spwvd'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    G  = tf_kernel_spwvd(N,param{k,:});
                    TF = filter_tfd(qtfd_wvd(z{k},N-1),G);
                    T(l,k) = toc;
                end
            end
        case 'ed'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    G  = tf_kernel_ed(N,param{k,:});
                    TF = filter_tfd(qtfd_wvd(z{k},N-1),G);
                    T(l,k) = toc;
                end
            end
        case 'bjd'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    G  = tf_kernel_bjd(N,param{k,:});
                    TF = filter_tfd(qtfd_wvd(z{k},N-1),G);
                    T(l,k) = toc;
                end
            end
        case 'bd'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    G  = tf_kernel_bd(N,param{k,:});
                    TF = filter_tfd(qtfd_wvd(z{k},N-1),G);
                    T(l,k) = toc;
                end
            end
        case 'mbd'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    G  = tf_kernel_mbd(N,param{k,:});
                    TF = filter_tfd(qtfd_wvd(z{k},N-1),G);
                    T(l,k) = toc;
                end
            end
        case 'embd'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    G  = tf_kernel_embd(N,param{k,:});
                    TF = filter_tfd(qtfd_wvd(z{k},N-1),G);
                    T(l,k) = toc;
                end
            end
        case 'ckd'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    G  = tf_kernel_ckd(N,param{k,:});
                    TF = filter_tfd(qtfd_wvd(z{k},N-1),G);
                    T(l,k) = toc;
                end
            end
        case 'rgd'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    G  = tf_kernel_rgd(N,tf2af(qtfd_wvd(z{k},N-1)),param{k,:});
                    TF = filter_tfd(qtfd_wvd(z{k},N-1),G);
                    T(l,k) = toc;
                end
            end
        case 'mdd'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    G  = tf_kernel_mdd(N,tf2af(qtfd_wvd(z{k},N-1)),param{k,:});
                    TF = filter_tfd(qtfd_wvd(z{k},N-1),G);
                    T(l,k) = toc;
                end
            end
        case 'dgf'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    TF = qtfd_dgf(qtfd_wvd(z{k},N-1),param{k,:});
                    T(l,k) = toc;
                end
            end
        case 'mpd'
            T = zeros(L,K);
            for k = 1:K
                for l = 1:L
                    tic;
                    Atoms = mp_atoms(z{k},dict,1);
                    TF = (N/size(Atoms,2)).*mpd(Atoms);
                    T(l,k) = toc;
                end
            end
    end
    save(['Data/Computation time/comtime_' idx '.mat'],'T');
end