function p2pModification(corpus, dtw, frameLength, overlap, savefolder, savefile, lpc)
% input
% corpus: using corpus
% maxRelErr: the target modifiction relative error
% savefolder: 
% savefile: 

    if nargin < 7, lpc = ''; end

    % selection of the measurement method for spectral tilt
    method = 6;
    % coefficient for pre/de emphasis filter
    emp_coeff = -0.97; 
    fs = 16000;
    frameSize = fix(frameLength * fs);
    win = hann(frameSize,"periodic");
    shiftSize = fix((1 - overlap) * frameSize);

    threshold = 0.1;
    
    result = table('Size', [height(corpus) * 250 7], ...
                    'VariableTypes', {'string','double','double','double','double','double','double'}, ...
                    'VariableNames',{'utterance','segmentIdx','tilt_p','tilt_l','tilt_m','relerr_frame','relerr_utter'});
    
    size = height(corpus);
    
    mkdir(savefolder);
    
    utter_start = 1;

    for i = 1 : size

        if mod(i,100) == 0, disp(i), end
        
        speaker = corpus.SPKR{i};
        utterance = corpus.UTTERANCE{i};
        fname_p = corpus.FNAME_P{i};
        fname_l = corpus.FNAME_L{i};
        
        if strcmpi(dtw,'lombard')
            audio_p = audioread(fullfile("lombardgrid_paired","audio", [fname_p, '.wav']));
            audio_l = audioread(fullfile("lombardgrid_paired","audio_dtw","lombard_alignedto_plain",...
                                [fname_l, '_aligned.wav']));
        elseif strcmpi(dtw,'plain')
            audio_p = audioread(fullfile("lombardgrid_paired","audio_dtw","plain_alignedto_lombard",...
                                [fname_p, '_aligned.wav']));
            audio_l = audioread(fullfile("lombardgrid_paired","audio", [fname_l, '.wav']));
        else
            disp('input plain or lombard to select using utterances modified by dtw');
            return;
        end

        [frames_p, tilt_p_current, voiced_p] = getFrameInfo(audio_p, fs, win, frameSize, shiftSize, threshold, method);
        [frames_l, tilt_l_current, voiced_l] = getFrameInfo(audio_l, fs, win, frameSize, shiftSize, threshold, method);

        voicedIdx = find(voiced_p & voiced_l);

        tilt_p = cell(length(voicedIdx),1);
        tilt_l = tilt_p;
        tilt_m = tilt_p;
        audio_m = zeros(length(audio_p),1);

        itr = length(voiced_p);
        k = 1;

        startFrameIdx = 1;
        endFrameIdx = startFrameIdx + frameSize - 1;

        for j = 1 : itr
            
            frame_p = frames_p(:,j);
            frame_l = frames_l(:,j);

            if k <= length(voicedIdx)

                if j == voicedIdx(k)
                
    
                    % Measure the spectral-tilts of plain/lombard segment
                    tilt_p{k} = tilt_p_current(voicedIdx(k));
                    tilt_l{k} = tilt_l_current(voicedIdx(k));
        
                    if strcmpi(lpc,'lpc')
                        
                        frame_m = LPCTransplantation(frame_p, frame_l, emp_coeff, fs);

                    else
        
                        frame_m = iterateSTfilter(frame_p, fs, tilt_p{k}, tilt_l{k}, emp_coeff, method);

                    end

%                     soundsc(frame_m,fs);
%                     pause(length(frame_m) * 1.5 / fs);
        
                    tilt_m{k} = myGetSpectralTilt(frame_m, fs, method);
                    k = k + 1;
    
                else
                    frame_m = frame_p;
                end
            end

            audio_m(startFrameIdx:endFrameIdx) = audio_m(startFrameIdx:endFrameIdx) + frame_m;
            % Update audio_m
     
            startFrameIdx = startFrameIdx + shiftSize;
            endFrameIdx = startFrameIdx + frameSize - 1;
        end

        % Gain match
        audio_m = audio_m * rms(audio_p)/rms(audio_m);
        soundsc(audio_m,fs);

        utter_end = utter_start + k - 2;
        
        result.utterance(utter_start : utter_end) = repmat({fname_p},length(voicedIdx),1);
        result.segmentIdx(utter_start : utter_end) = 1:length(voicedIdx);;
        result(utter_start : utter_end, {'tilt_p','tilt_l','tilt_m'}) = [tilt_p tilt_l tilt_m];

        tilt_err_frame = (cell2mat(tilt_l) - cell2mat(tilt_m)) ./  cell2mat(tilt_l);
        result{utter_start:utter_end, 'relerr_frame'} = tilt_err_frame;
        
        result{utter_start:utter_end, 'relerr_utter'} = mean(abs(tilt_err_frame));
    
        % Store modified utterance
        audioName = fullfile(savefolder, [speaker '_m_' utterance '.wav']);
        audiowrite(audioName, audio_m, fs);

        utter_start = utter_end + 1;
 
    end
       
    % remove empty rows
    result = rmmissing(result);
    save(fullfile(savefolder,[savefile, '.mat']),'result','-mat');

end

function [frames, tilt, voiced] = getFrameInfo(audio, fs, win, frameSize, shiftSize, threshold, method)

        weight = 0.8;
        
        itr = fix((length(audio) - (frameSize -  shiftSize)) / shiftSize) - 1;
        
        startIdx = 1;
        endIdx = startIdx + frameSize - 1;
        
        tilt = zeros(itr,1);
        intensity = zeros(itr,1);

        frames = zeros(frameSize, itr);
        
        for i = 1:itr
        
            frames(:,i) = audio(startIdx : endIdx) .* win;
            frames(:,i) = frames(:,i) - mean(frames(:,i));
        
            tilt(i) = myGetSpectralTilt(frames(:,i), fs, method);
            intensity(i) = rms(frames(:,i));
        
            startIdx = startIdx + shiftSize;
            endIdx = startIdx + frameSize - 1;
        
        end
        
        %%Normalize
    
        spectilt = normalize(tilt, 'range');
        intensity = normalize(intensity, 'range');
         
        detection = intensity - weight * spectilt;
        voiced = detection > threshold;

        voicedforplot = resample(double(voiced), length(audio), itr);

        figure;
        plot(audio);
        hold on
        plot(voicedforplot);
        hold off
end