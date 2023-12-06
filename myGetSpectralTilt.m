% Script for the computation of spectral tilt.
%
% Arguments:
%   x: a vector of samples
%   sr: sampling rate
%   method:
%    1: as the difference between two adjacent spectral bands (40~1k : 1k~12.5k)
%    2: as the lin-reg of the PSD of x, resampled @ 16000
%    3: as the lin-reg of the LPC
%    4: as the lin-reg of the LPC residual
%    5: as the lin-reg of the magnitude spectrum of x
%    6: as the lin-reg of 1/3 energy bands (by Cassia Valentini)
%   preemph: using preemphasis (1) or not (0)
%   dbX8ve: return results in dB per Octaves (1) or per spectrum (0)
%
% Returns an array 1x2 [slope y0]
%
% Created by Julian Villegas
% Lista project 2011
%%
function result=myGetSpectralTilt(x,sr,method,emp_coeff)

if nargin<4
    emp_coeff = nan;
end

if ~isnan(emp_coeff) 
    x=filter([1 emp_coeff],1,x);
end

%%If there is not any signal in x, tilt=0;
if db(rms(x)) == -Inf
    result = 0;
    return;
end

% Remove DC  offset
x = x - mean(x);

switch method
    case 1
%         result = spectralSlope(x,sr,'Window',round(0.02*sr),'OverlapLength');
        
    case 2
        % as the lin-reg of the PSD of x, resampled @ 16000
        len = 1024;                     % FFT length
        m = len/2;                      % number of distinct frequency bins
        dt = 0.001;                     % size of time step
        fc = 1/(2 * dt);                % the critical frequency
        fb = fc * (0:m)/m;              % the freq bins (513)
        X = fft(x, len);                % FFT (contains neg/pos freq) 
        Pxx =  log10(X .* conj(X) / len);       % PSD = |Y|^2
        Pxx(m+2:len) = [];
%         Pxx(2:m+1) = 2 * Pxx(2:m+1);    % compensate for missing neg freq
        result = polyfit(fb', Pxx, 1);  
        result = 1000 * result(1);


    case {3,4}
        % as the lin-reg of the LPC
        lpcOrder = 20;
%         freqs = 0:0.01:(sr/2);
        freqs = 500:0.01:7500;
        ax = real(lpc(x,lpcOrder));
        est = freqz(1,ax,freqs,sr);     %estimated spectrum
        axSpec = 20*log10(abs(est));    %calc magnitude spectrum
        result = polyfit(freqs, axSpec, 1);

        if method == 4
            glottalX = filter([1 ax(2:end)],1,x);
            ax = real(lpc(glottalX,lpcOrder));
            axSpecGlot = 20*log10(abs(freqz(1,ax,freqs,sr)));
            result = polyfit(freqs,axSpecGlot,1);
        end

        result = 1000*result(1); 
    
    case 5
        % as the lin-reg of the magnitude spectrum
        result = 0;

    case 6
        % as the lin-reg of 1/3 energy bands (by Cassia Valentini)
        
        x = resample(x,16000,sr);
        N_window = round(.02*sr); % 20 msec window (note default overlap in spectrogram is 50%)
        [Ls] = ltas(x,N_window,sr,1);
        nfft = 2*(length(Ls)-1);
%         fft_bins = [0:(.5*nfft)] * sr/nfft;
        
        % Center frequencies according to one third octave band  
        bcf = [500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000];
%         bcf = [160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000];
        Nb = length(bcf);
        
        fl = [450 565 715 900 1125 1425 1800 2250 2825 3575 4500 5650 7150];
%         fl = [100 180 225 283 358 450 565 715 900 1125 1425 1800 2250 2825 3575 4500 5650 7150];
        fu = [fl(2:end) 8000];
        
        % finding the indices corresponding the upper and lower band frequencies
        fli=max(1,round(nfft*fl/sr));
        fui=round(nfft*fu/sr);
        fui(end) = nfft/2;
        
        E=zeros(1,Nb);
        for n=1:length(fli)
            E(n) = sum(10.^(Ls(fli(n):fui(n))/10)); % estimating the spectrum level in each of the individual bands
            E(n) = E(n)/(fui(n)-fli(n)+1);
            E(n) = 10*log10(E(n));
        end
        
        % Linear regression
        result = polyfit(1:Nb,E,1);
    
        stfit = polyval(result,1:Nb);
        result = (stfit(end)-stfit(1))/4; % dB/octave
%         result = (stfit(18)-stfit(3))/5; % dB/octave
        

    case 7
        % as the c1 (cepstrum coefficient)
        n = 2^nextpow2(length(x));

        dft = abs(fft(x,n));

        cps = -ifft(20*log10(dft));

        result = cps(2);

end



function [Ls,w,Ps] = ltas(s,N,Fs,type)
% Calculates the Long term average spectrum
%
% Inputs:
%     s - signal (col vector)
%     N - window size for fft
%     Fs - sampling rate (Hz)
%     type - 0: dB units / 1: SPL dB units
% Outputs:
%     Ls - (col vector)
%     w  - frequency values (col vector)
% Author: Cassia Valentini-Botinhao
% Created: 06.06.10

    [S,w,T,Ps] = spectrogram(s,N,[],[],Fs);
    Ps = mean(Ps,2);
    
    Po = 20*10^-6; % reference pressure
    
    if type == 1
        Ps = Ps./(Po.^2);
    end
    
    Ls = 10*log10(Ps);
    end

end