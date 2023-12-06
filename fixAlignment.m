%% 
% ・Exclude voiceless segments: We are only interested in the modification for 
% voiced sounds
% 
% ・Detect the beginning of utterances: Grid corpus has the alignements, but 
% we need to give them offset
% 
% ・To detect first voiced sound, we use spectral tilt and intensity
% 
% 

clearvars
load lombardgrid_paired/corpusCleaned.mat

[~, fs] = audioread(fullfile('lombardgrid_paired/audio/',...
                [corpusCleaned.FNAME_P{1}, '.wav']));

frameSize = fix(0.02 * fs);
win = hann(frameSize,"periodic");
overlap = 0.50;
shiftSize = fix((1 - overlap) * frameSize);
threshold_tilt = 0.1;
threhsold_intensity = 0.1

%%  Assume that each utterances has sounds less than 20 
size = [height(corpusCleaned) * 30*2 4];

alignments = table('Size', size, 'VariableTypes', {'string','string','int32','int32'}, ...
             'VariableNames', {'utter_info','phone','duration','offset'});
%%
tic

endIdx = 0;

for i = 1:height(corpusCleaned)

    if ~mod(i,100), disp(i), end

    startIdx = endIdx + 1;


    [audio_p, fs] = audioread(fullfile('lombardgrid_paired/audio/',...
              [corpusCleaned.FNAME_P{i}, '.wav']));
    [audio_l, ~] = audioread(fullfile('lombardgrid_paired/audio/',...
              [corpusCleaned.FNAME_L{i}, '.wav']));

    alignment_p = getAlignment(corpusCleaned.FNAME_P{i},fs);
    alignment_l = getAlignment(corpusCleaned.FNAME_L{i},fs);

    alignment_p_fixed = alignmentOffset(audio_p, alignment_p, fs, frameSize, win, overlap, threshold);
    alignment_l_fixed = alignmentOffset(audio_l, alignment_l, fs, frameSize, win, overlap, threshold);

    endIdx = startIdx + height(alignment_p_fixed) - 1;
    alignments(startIdx : endIdx, :) = alignment_p_fixed;

    startIdx = endIdx + 1;
    endIdx = startIdx + height(alignment_l_fixed) - 1;
    alignments(startIdx : endIdx, :) = alignment_l_fixed;

    if alignment_p_fixed.offset(end)+alignment_p_fixed.duration(end)  > length(audio_p)
        disp(alignment_p_fixed.utter_info(1));
    end

    if alignment_l_fixed.offset(end)+alignment_l_fixed.duration(end)  > length(audio_l)
        disp(alignment_l_fixed.utter_info(1));
    end
    

end

toc

alignments = rmmissing(alignments);

save(fullfile('lombardgrid_paired','alignments_fixed'),"alignments");
%%
function alignment = alignmentOffset(audio, alignment, fs, frameSize, win, overlap, threshold)
    
        method = 6;
        emp_coeff = -0.97;
        shiftSize = fix((1 - overlap) * frameSize);
        voiceless = {'ch_','f_','k_','p_','s_','sh_','t_','th_'};
        
        itr = fix((length(audio) - (frameSize -  shiftSize)) / shiftSize) - 1;
        
        startIdx = 1;
        endIdx = startIdx + frameSize - 1;
        
        spectilt = zeros(itr,1);
        intensity = zeros(itr,1);
        
        for i = 1:itr
        
            audioFrame = audio(startIdx : endIdx) .* win;
            audioFrame = audioFrame - mean(audioFrame);
        
            spectilt(i) = myGetSpectralTilt(audioFrame, fs, method);
            intensity(i) = rms(audioFrame);
        
            startIdx = startIdx + shiftSize;
            endIdx = startIdx + frameSize - 1;
        
        end
        
        %%Normalize
    
        spectilt = normalize(spectilt, 'range');
        intensity = normalize(intensity, 'range');
        
        spectilt = resample(spectilt, length(audio), itr);
        intensity = resample(intensity, length(audio), itr);
    
    %% Intensity is more strong parameter of voiced part than spectral tilt 
    %% Therefore, arbitrary weighting is done to detect voiced part
    % in the case of [l], spectilt should have smaller weighting than others
    % because [l] is not fricative/bilabial sound.
    if contains(alignment.phone{1},'l')
        detection = intensity - 0.15 * spectilt;
    elseif contains(alignment.phone{1},'b')
        detection = intensity - 0.5 * spectilt;
    elseif contains(alignment.phone{1}, 'p')
        detection = intensity - 0.6 * spectilt;
    else
        detection = intensity - 0.6 * spectilt;
    end
    
    %% In the voiced sound, intensity > spectilt (both normalized)
    firstZero = find(detection > threshold, 1) + frameSize/2;
    firstVoiced = find(~contains(alignment.phone, voiceless),1);
    
    alignment.offset = alignment.offset - alignment.offset(firstVoiced) + firstZero;


end

function alignment = getAlignment(fname, fs)

%     voiceless = {'SIL','ch_','f_','k_','p_','s_','sh_','t_','th_'};
    
    fileID = fopen(fullfile('lombardgrid_paired/alignment/',...
                  [fname, '.json']));
    raw = fread(fileID, inf);
    str = char(raw');
    fclose(fileID);
    
    alignment = struct2table(jsondecode(str).(fname));
    % Exclude silence and voiceless sounds from the alignment
    alignment = alignment(~contains(alignment.phone,'SIL'),:);
    
    utter_info = table(repmat({fname}, height(alignment), 1) ...
            ,'VariableNames', {'utter_info'});
    
    alignment = [utter_info alignment(:,{'phone','duration','offset'})];

    % To ease the process after this, convert from seconds to indices
    alignment{:,{'duration','offset'}} = floor(alignment{:,{'duration','offset'}} * fs);

end

function [] = plotAlignment(audio, alignment, fs)
    
    t = linspace(1,length(audio)/fs, length(audio));

    figure;
    plot(t, audio);
    hold on
    xline(t(alignment.offset), 'r:', 'LineWidth', 1.5);
    xline(t(alignment.offset + alignment.duration), 'r-.', 'LineWidth',1.5);

    xline(t(alignment.offset(1)), 'r', 'LineWidth', 2);
    xline(t(alignment.offset(end) + alignment.duration(end)), 'r', 'LineWidth', 2);
    
    title(sprintf('Result of fixAlignment: %s',alignment.utter_info{1}));

    hold off

end