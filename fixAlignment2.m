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
load lombardgrid_paired/vadThreshold.mat

% tilt_threshold = tilt_threshold - abs(tilt_threshold * 0.2);
% intensity_threshold = intensity_threshold + abs(intensity_threshold * 0.);

[~, fs] = audioread(fullfile('lombardgrid_paired/audio/',...
                [corpusCleaned.FNAME_P{1}, '.wav']));

shiftSize = fix(0.02 * fs);

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
    
    spkr = corpusCleaned.SPKR{i};

    [alignment_p_fixed, alignment_l_fixed] ...
    = alignmentOffset(audio_p, audio_l, alignment_p, alignment_l, fs, shiftSize);


    endIdx = startIdx + height(alignment_p_fixed) - 1;
    alignments(startIdx : endIdx, :) = alignment_p_fixed;

    startIdx = endIdx + 1;
    endIdx = startIdx + height(alignment_l_fixed) - 1;
    alignments(startIdx : endIdx, :) = alignment_l_fixed;

end

toc

alignments = rmmissing(alignments);

save(fullfile('lombardgrid_paired','alignments_fixed'),"alignments");


%%
function [alignment_p, alignment_l] = alignmentOffset(audio_p, audio_l, alignment_p, alignment_l, fs, shiftSize)
   
    alignment_p.offset = alignment_p.offset - alignment_p.offset(1) + 1;
    alignment_l.offset = alignment_l.offset - alignment_l.offset(1) + 1;

    itr_p = round((length(audio_p) - alignment_p.duration(end) - alignment_p.offset(end)) / shiftSize * 0.9);
    itr_l = round((length(audio_l) - alignment_l.duration(end) - alignment_l.offset(end)) / shiftSize * 0.9);

    alignment_p_tmp = alignment_p;
    alignment_l_tmp = alignment_l;

    startIdx = [-Inf, nan, nan];

    for i = 1:itr_p
        alignment_p_tmp.offset = alignment_p.offset + shiftSize * (i-1);

        for j = 1:itr_l
            alignment_l_tmp.offset = alignment_l.offset + shiftSize * (j-1);
            [zerox_score, rms_score] = measureVoiced(audio_p, audio_l, alignment_p_tmp, alignment_l_tmp, fs);

            if rms_score * zerox_score > startIdx(1)
                startIdx = [zerox_score * rms_score, i, j];
            end

        end



    end
    
%     alignment_p.offset = alignment_p.offset + shiftSize * relerr_min(2);
%     alignment_l.offset = alignment_l.offset + shiftSize * relerr_min(3);
    alignment_p.offset = alignment_p.offset + shiftSize * startIdx(2);
    alignment_l.offset = alignment_l.offset + shiftSize * startIdx(3);

end






function [zerox_score, rms_score] = measureVoiced(audio_p, audio_l, alignment_p, alignment_l, fs)

    % selection of the measurement method for spectral tilt
    method = 6;
    emp_coeff = -0.97;
    voiceless = {'ch_','f_','k_','p_','s_','sh_','t_','th_'}; 

    plainLength = height(alignment_p);
    
    %% Rounding of segment ends by 30ms by hann window
    edgeFrame = fix(0.03 * fs);

    rms_score = 0;
    zerox_score = 0;

    segStart_p = alignment_p.offset + edgeFrame;
    segEnd_p = alignment_p.offset + alignment_p.duration + edgeFrame - 1;
    segStart_l = alignment_l.offset + edgeFrame;
    segEnd_l = alignment_l.offset + alignment_l.duration + edgeFrame - 1;

    % Add 30ms silence to the edges of the audio signal to do windowing 
    audio_p = [zeros(edgeFrame, 1); audio_p; zeros(edgeFrame, 1)];
    audio_l = [zeros(edgeFrame, 1); audio_l; zeros(edgeFrame, 1)];
                    
    for i = 1 : plainLength

        seg_p = audio_p( (segStart_p(i) - edgeFrame * 0.5) : (segEnd_p(i) + edgeFrame * 0.5) );  
        seg_l = audio_l( (segStart_l(i) - edgeFrame * 0.5) : (segEnd_l(i) + edgeFrame * 0.5) );

        seg_p = seg_p - mean(seg_p);
        seg_l = seg_l - mean(seg_l);

        seg_p = edgeWindowing(seg_p, edgeFrame * 2);
        seg_l = edgeWindowing(seg_l, edgeFrame * 2);

        seg_p = filter([1 emp_coeff], 1, seg_p);
        seg_l = filter([1 emp_coeff], 1, seg_l);

        if ~contains(alignment_p.phone{i},voiceless)
            zerox_tmp = (1 - zerocrossrate(seg_p)) + (1 - zerocrossrate(seg_l));
            rms_tmp = rms(seg_p) + rms(seg_l);

            zerox_score = zerox_score + zerox_tmp;
            rms_score = rms_score + rms_tmp;
        else
            zerox_tmp = zerocrossrate(seg_p) + zerocrossrate(seg_l);
            rms_tmp = rms(seg_p) + rms(seg_l);

            zerox_score = zerox_score + zerox_tmp * 0.5;
            rms_score = rms_score + rms_tmp * 0.5;

        end
    end
        
    rms_tmp = rms(audio_p(1:segStart_p(1))) + rms(audio_l(1:segStart_l(1)))+...
              rms(audio_p(segEnd_p(end) : end)) + rms(audio_l(segEnd_l(end) : end));

    rms_score = rms_score - rms_tmp;


end







function alignment = getAlignment(fname, fs)

%     voiceless = {'SIL','ch_','f_','k_','p_','s_','sh_','t_','th_'};
    
    fileID = fopen(fullfile('lombardgrid_paired/alignment/',...
                  [fname, '.json']));
    raw = fread(fileID, inf);
    str = char(raw');
    fclose(fileID);
    
    alignment = struct2table(jsondecode(str).(fname));
    % Exclude silence segment from the alignment
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