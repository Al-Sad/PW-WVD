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
%                Piece-Wise Spline Wigner-Ville Distribution
%
%  Syntax : [Lz, Lp] = pw_wvd(IF,IP,IA,W,fs)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% IF : The signal IF (P x N).
% IP : The signal IP (P x N).
% IA : The signal IA (P x N).
% W  : The signal time-support (P x N).
% fs : Sampling frequency in Hz.
%
% <OUTPUTs>
% Lz : PW-WVD of the signal (N x N).
% Lp : PW-WVD of each signal component.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [Lz, Lp] = pw_wvd(IF,IP,IA,W,fs)
%% Initialization
[P, N] = size(IF);
Lz     = zeros(N,N);
Lp     = cell(1,P);

%% Parameters
t   = 0:1/fs:(N-1)/fs; % time array in seocnds
Q   = sum(W,2)-1;
nn  = 1:N;
W   = logical(W);
K   = 2*N;
L   = (N-1)*~mod(N,2) + (N)*mod(N,2);
Lh  = fix(L/2);
tau = -Lh:Lh;
mm  = 1 + rem(N + tau, N);

%% Main
for p = 1:P
    % Remove the IF samples that are outside the time support
    IF(p,~W(p,:)) = nan;
    % Approximate the pth IF using linear B-spline
    n = nn(W(p,:));
    dt = diff(t(W(p,:)));
    df = diff(IF(p,W(p,:)));
    m  = df./dt;
    b  = IF(p,n(1:end-1)) - m.*t(n(1:end-1));
    pw_if = (t'*m + repmat(b,N,1))'.*W(p,:);
    % Approximate the pth IP and fix its boundary conditions
    pw_ip = (0.5*(t.^2)'*m + t'*b)';
    pw_ip(1,:) = pw_ip(1,:) - (pw_ip(1,n(1)) - IP(p,n(1)));
    for q = 2:Q(p)
        c = (pw_ip(q,n(q)) - pw_ip(q-1,n(q)));
        pw_ip(q,:) = pw_ip(q,:) - c;
    end
    % Compute the analytic signal using the piece-wise LFM model
    s      = IA(p,:).*exp(1j*2*pi*pw_ip).*(pw_if > 0 & pw_if < fs/2).*W(p,:);
    z      = fft(s,N,2);           % Transforming S1 into the frequency domain
    z(:,1) = z(:,1);               % Holding the zero frequency amplitude/energy
    z(:,2:N/2) = 2*z(:,2:N/2);     % Accepting and doubling half of the frequencies
    z(:,N/2+1:N) = 0;              % Rejecting half of the frequencies
    z = [ifft(z,[],2) zeros(Q(p),N)]; % Padding with zeros
    % Compute the WVD for each q spline
    KTL  = zeros(K,K);
    for i = 1:n(1)
        Z1  = z(1,1 + rem(K + i - 1 + tau, K));
        Z2  = z(1,1 + rem(K + i - 1 - tau, K));
        KTL(i,mm) = Z1.*conj(Z2);
    end
    for q = 2:Q(p)
        Z1  = z(q,1 + rem(K + n(q) - 1 + tau, K));
        Z2  = z(q,1 + rem(K + n(q) - 1 - tau, K));
        KTL(n(q),mm) = Z1.*conj(Z2);
    end
    for i = n(Q(p)+1):N
        Z1  = z(1,1 + rem(K + i - 1 + tau, K));
        Z2  = z(1,1 + rem(K + i - 1 - tau, K));
        KTL(i,mm) = Z1.*conj(Z2);
    end
    Wz = real(fft(KTL,N,2));
    Lz(1:N-1,:) = Wz(1:N-1,:) + Lz(1:N-1,:);
    Lp{p} = Wz(1:N-1,:);
end
end