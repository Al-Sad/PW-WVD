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
% This main script generates the database signals and parameters and saves
% them in pw_wvd_database.mat. It produces 1000 multi-component
% non-stationary signals sampled at 1 Hz and characterized by random number
% of components between 1 and 4, polynomial IF laws with random orders
% between 1 and 3, random constant IA between 0.5 and 1, and random time
% support within 0 and 255. The signal random time support is constrained
% with minimum and maximum IF curve lengths to produce signals with
% realistic durations.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
fs   = 1;
N    = 256;
M    = 256;
Pmax = 4;
Rmax = 3;
K    = 1000;

%% Generate P and R for all K signals
F = [];
for p = 1:Pmax
    F = [F ; zeros(nchoosek(p+Rmax-1,p),Pmax-p) combs_rep(Rmax,p)];
end
Ls = size(F,1);
S = repmat(1:Ls,1,floor(K/Ls));
S = [S randi(size(F,1),1,K- Ls*floor(K/Ls))];
S = S(randperm(length(S)));
PR = F(S,:);

%% Signal Generation and PW-WVD computation
z   = cell(K,1);
Lz  = cell(K,1);
IF  = cell(K,1);
IP  = cell(K,1);
IA  = cell(K,1);
Wp  = cell(K,1);
rng(1);
for k = 1:K
    R = PR(k,PR(k,:)~=0);
    [z{k}, IF{k}, IP{k}, IA{k}, Wp{k}] = signal_generator(R,fs,N,M);
    disp(100*k/K)
end

%% Saving
save(['Data' filesep 'Database' filesep 'pw_wvd_database.mat'],'z','IF','IP','IA','Wp','PR','S','fs','N','K');
