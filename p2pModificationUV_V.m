function p2pModificationUV_V(corpus, alignments, frameLength, overlap, savefolder, savefile, method, lpc)
% input
% corpus: using corpus
% maxRelErr: the target modifiction relative error
% savefolder: 
% savefile: 

    if nargin < 8, lpc = ''; end

    % selection of the measurement method for spectral tilt
    
    result = table('Size', [height(corpus) * 400 7], ...
                    'VariableTypes', {'string','double','double','double','double','double','double'}, ...
                    'VariableNames',{'utterance','segmentIdx','tilt_p','tilt_l','tilt_m','abserr_frame','abserr_utter'});
    
    size = height(corpus);
    
    mkdir(savefolder);
    
    utter_start = 1;
    fs = 16000;

    for i = 1 : size

        if mod(i,100) == 0, disp(i), end
        
        speaker = corpus.SPKR{i};
        utterance = corpus.UTTERANCE{i};
        fname_p = corpus.FNAME_P{i};
        fname_l = corpus.FNAME_L{i};
        
        
        audio_p = audioread(fullfile("lombardgrid_paired","audio", [fname_p, '.wav']));
        audio_l = audioread(fullfile("lombardgrid_paired","audio_dtw","lombard_alignedto_plain",...
                                [fname_l, '_aligned.wav']));

        alignment = alignments(strcmp(alignments.utter_info, fname_p),:);

        [audio_m, tilt_p, tilt_l, tilt_m, abserr_frame, abserr_utter] = modify(audio_p, audio_l, alignment, frameLength, overlap, lpc, fs, method);

        itr = length(tilt_p);

        utter_end = utter_start + itr - 1;

            
        result.utterance(utter_start : utter_end) = repmat({fname_p},length(itr),1);
        result.segmentIdx(utter_start : utter_end) = 1:itr;;
        result(utter_start : utter_end, {'tilt_p','tilt_l','tilt_m'}) = [tilt_p tilt_l tilt_m];
        result{utter_start:utter_end, 'abserr_frame'} = abserr_frame;
        result{utter_start:utter_end, 'abserr_utter'} = abserr_utter;
    
        % Store modified utterance
        audioName = fullfile(savefolder, [speaker '_m_' utterance '.wav']);
        audiowrite(audioName, audio_m, fs);

        utter_start = utter_end + 1;
 
    end
       
    % remove empty rows
    result = rmmissing(result);
    save(fullfile(savefolder,[savefile, '.mat']),'result','-mat');

end


function [audio_m, tilt_p, tilt_l, tilt_m, abserr_frame, abserr_utter] = modify(audio_p, audio_l, alignment, frameLength, overlap, lpc, fs, method)

        % coefficient for pre/de emphasis filter
        emp_coeff = -0.97; 
        frameSize = round(frameLength * fs);
        win = hann(frameSize,"periodic");
        shiftSize = round((1 - overlap) * frameSize);

        % Modify utterance segment only
        utterStart = alignment.offset(1);
        utterEnd = alignment.offset(end) + alignment.duration(end) - 1;

        audio_m = audio_p;


        itr = round(((utterEnd - utterStart)- (frameSize -  shiftSize)) / shiftSize) - 1;

        tilt_p = cell(itr,1);
        tilt_l = tilt_p;
        tilt_m = tilt_p;

        startFrameIdx = utterStart;
        endFrameIdx = startFrameIdx + frameSize - 1;

        for j = 1 : itr

            frame_p = audio_p(startFrameIdx:endFrameIdx) .* win;
            frame_l = audio_l(startFrameIdx:endFrameIdx) .* win;
            
            % Measure the spectral-tilts of plain/lombard segment
            tilt_p{j} = myGetSpectralTilt(frame_p, fs, method);
            tilt_l{j} = myGetSpectralTilt(frame_l, fs, method);


            if db(rms(frame_p)) > -Inf && db(rms(frame_l)) > -Inf
    
                if strcmpi(lpc,'lpc')
                    
                    frame_m = LPCTransplantation(frame_p, frame_l, emp_coeff, fs);
    
                else
    
                    frame_m = iterateSTfilter(frame_p, fs, tilt_p{j}, tilt_l{j}, method);
    
                %audio_voiced(startFrameIdx:endFrameIdx) = audio_voiced(startFrameIdx:endFrameIdx) + frame_m;
    
                end
    
                tilt_m{j} = myGetSpectralTilt(frame_m, fs, method);

            else

                frame_m = frame_p;
                tilt_m{j} = tilt_l{j};
            
            end

            audio_m(startFrameIdx:endFrameIdx) = audio_m(startFrameIdx:endFrameIdx) .* (1-win) + frame_m;

            % Update audio_m
     
            startFrameIdx = startFrameIdx + shiftSize;
            endFrameIdx = startFrameIdx + frameSize - 1;
        end

        % Gain match
        audio_m = audio_m * rms(audio_p)/rms(audio_m);


        abserr_frame = cell2mat(tilt_l) - cell2mat(tilt_m);
        abserr_utter = mean(abs(abserr_frame));

end





function [frames, tilt, voiced] = getFrameInfo(audio, fs, win, frameSize, shiftSize, threshold, method)

        weight = 0.4;
        
        itr = round((length(audio) - (frameSize -  shiftSize)) / shiftSize) - 1;
        
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

%         voicedforplot = resample(double(voiced), length(audio), itr);
% 
%         figure;
%         plot(audio);
%         hold on
%         plot(voicedforplot);
%         hold off
end