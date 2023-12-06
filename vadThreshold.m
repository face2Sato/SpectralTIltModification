clearvars;
load lombardgrid_paired/corpusCleaned.mat
% load utterance_tilt_intensity.mat
load utterance_tilt_intensity.mat

tilt_threshold = zeros(55,2);
intensity_threshold = zeros(55,2);



targetTilt = -6;
targetIntensity = -60;

for i = 2:55
    
    spkr_p = utterance_info(contains(utterance_info.utterance, sprintf('s%d_p_s',i)),:);
    spkr_l = utterance_info(contains(utterance_info.utterance, sprintf('s%d_l_s',i)),:);
    
    [tilt_threshold(i,1), intensity_threshold(i,1)] = getpeak(spkr_p, targetTilt, targetIntensity);

    [tilt_threshold(i,2), intensity_threshold(i,2)] = getpeak(spkr_l, targetTilt, targetIntensity);    

end

save('lombardgrid_paired/vadThreshold',"tilt_threshold","intensity_threshold");



function [tilt_threshold, intensity_threshold] = getpeak(spkr, targetTilt, targetIntensity)

    intensity = spkr.intensity(-80 < spkr.intensity & spkr.intensity < -40);
    tilt = spkr.tilt(-12 < spkr.tilt & spkr.tilt < 2);

    fs = 1024;
    
    windowSize = fs /8; 
    b = (1 / windowSize) * ones(1, windowSize);
    a = 1;
    
    [value, bins] = histcounts(tilt,fs);
    bins = (bins(1:end-1) + bins(2:end)) / 2;
    value = filter(b, a, value);
    % value = h1.Values;
    [~,locs] = findpeaks(-value,"MinPeakDistance",fs/3.5);

    [~, Idx] = min(abs(bins(locs) - targetTilt));

    tilt_threshold = bins(locs(Idx));

    figure;
    plot(bins, value);
    xline(bins(locs(Idx)));


    [value, bins] = histcounts(intensity,fs);
    bins = (bins(1:end-1) + bins(2:end)) / 2;
    value = filter(b, a, value);
    [~,locs] = findpeaks(-value,"MinPeakDistance",fs/3.5);

    [~, Idx] = min(abs(bins(locs) - targetIntensity));

    figure;
    plot(bins, value);
    xline(bins(locs(Idx)));


    intensity_threshold = bins(locs(Idx));

end
