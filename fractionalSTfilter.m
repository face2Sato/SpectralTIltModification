function y= fractionalSTfilter(x,sr,alpha,f0,r,N)
% modifyTilt filters a signal with a filter that has a specified spectral tilt
%segment
% SYNOPSIS: y=modifyTilt(x,sr,alpha,f0,f1,r,N)
%
% INPUT x: monophonic audio
%       sr: sampling rate
%       alpha: (opt) spectral tilt. Default: -0.5 (pink noise)
%       f0: (opt) lowest frequency. Default: 20 Hz
%       r:  (opt) interval between poles and zeros. Default: 6/5 (minor third)
%       N:  (opt) number of poles. Default: 39
%
% OUTPUT y: the filtered signal
%
% REMARKS This code is based on the code presented in "Closed Form Fractional
% Integration and Differentiation via Real Exponentially Spaced Pole-Zero Pairs"
% by Julius Orion Smith and Harrison Freeman Smith and the implementation
% found in Faust Programming language
%
% SEE ALSO testHelpFunction
%
% AUTHOR    : Julian Villegas
% $DATE     : 02-Jun-2021 09:13:56 $
% $Revision : 1.00 $
% DEVELOPED : 9.10.0.1669831 (R2021a) Update 2
% FILENAME  : hoge.m

if nargin<6
    N = 40-1; % Number of poles: 4 m3 in an octave, 10 octaves minus one.
end
if nargin<5
    r = 6/5; % minor thirds
end
if nargin<4
    f0 = 20;
end
if nargin<3
    alpha = -1;
end
if nargin<2
    sr = 48000;
end
if nargin<1
    x = rand(sr,1);
    warning('insuficient arguments')
end

K = 3;
% The number of adding pole/zero pairs outside the band of interest

sr_tmp = 44000;

% Upsampling original sampling rate to 44000 to add pole-zero pairs outside the band of interest
% to approximate desired slope

x_len = length(x);

x = resample(x, sr_tmp, sr);

pzPlacement = -K : (N-1 + K);

w0 = 2 * pi * f0; % angular frequency

mp = w0*r .^ pzPlacement; % minus the poles
mz = w0*r .^ (-alpha + pzPlacement); % minus the zeros


% prewarping for bilinear transform
mph = prewarp(mp,sr_tmp,w0);
mzh = prewarp(mz,sr_tmp,w0);

for i=1:(N + 2 * K)
%     g = mph(i)/mzh(i);    
    [b1, b0, a1] = tf1s(1.0,mzh(i),mph(i),1,sr_tmp);
    x = filter([b0, b1],[1 a1], x);
end

%calculate the unity gain
g = mph(end) / mzh(end);
x = g .* x;

x = resample(x, sr, sr_tmp);

x = x(1:x_len);

y = x-mean(x);

%getTilt(y,sr)
end

function warped = prewarp(w,sr,wp)
    T = 1/sr;
    warped = wp * tan(w*T/2)/tan(wp*T/2);
end

% tf1s(b1,b0,a0,1)

%        b1 s + b0
% H(s) = ----------
%           s + a0

function [b1d, b0d, a1d] = tf1s(b1,b0,a0,w1,sr)
c   = 1/tan(w1*0.5/sr); % bilinear-transform scale-factor, w1=1
d   = a0 + c;
b1d = (b0 - b1*c) / d;
b0d = (b0 + b1*c) / d;
a1d = (a0 - c) / d;
end


%% memo
% the amplitude of the output seems not to be normalized
% When alpha = 1, the amplitude of the output is quite big
% In comparison, when alpha = -1, the amplitude of the output is too small
% The function may not be work where about sr < 41000... I don't know why;;

