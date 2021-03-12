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
%                      Chirplet dictionary parameters
%
%  Syntax : parameters = chirplet_dictionary(N,fs)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% N  : The signal total number of time samples.
% fs : Sampling frequency in Hz.
%
% <OUTPUTs>
% parameters : The chirplet atoms parameters (4 x 5767168).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function parameters = chirplet_parameters(N,fs)
%% Adjustable dictionary size parameters
Ns = 11;
Nf = 64;
Nt = 128;
Nb = 64;

%% Atom parameters
scales  = 2.^linspace(10,0,Ns);
freqs   = linspace(0,fs,Nf).*pi;
tShifts = linspace(0,(N-1)/fs,Nt);
beta    = linspace(-fs,fs,Nb).*(pi/50);

%% Initialize the dictionary parameters
parameters = zeros(4,length(scales)*length(freqs)*length(tShifts)*length(beta));

%% Compute the dictionary parameters
cnt = 1;
for tau = 1:length(tShifts)
    for f = 1:length(freqs)
        for s = 1:length(scales)
            for b = 1:length(beta)
                parameters(:,cnt) = [scales(s) tShifts(tau) freqs(f) beta(b)];
                cnt = cnt + 1;
            end
        end
    end
end
end
