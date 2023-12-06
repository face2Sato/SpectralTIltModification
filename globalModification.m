function globalModification(corpus, alignments, savefolder, savefile, method, lpc)
% input
% corpus: using corpus
% maxRelErr: the target modifiction relative error
% savefolder: 
% savefile: 

    if nargin < 6, lpc = ''; end

    % coefficient for pre/de emphasis filter
    emp_coeff = -0.97; 
    fs = 16000;
    
    result = table('Size', [height(corpus) 5], ...
                    'VariableTypes', {'string','double','double','double','double'}, ...
                    'VariableNames',{'utterance','tilt_p','tilt_l','tilt_m','abserr_utter'});
    
    error = table('Size', [height(corpus) 2], 'VariableTypes',{'string','double'}, ...
                                               'VariableNames',{'utterance','error_type'});

    size = height(corpus);
    mkdir(savefolder);
    errorIndex = 1;

    %% Rounding of segment ends by 30ms by hann window
    edgeFrame = fix(0.03 * fs);

    for i = 1 : size
    
        if mod(i,100) == 0, disp(i), end
        
        speaker = corpus.SPKR{i};
        utterance = corpus.UTTERANCE{i};
        fname_p = corpus.FNAME_P{i};
        fname_l = corpus.FNAME_L{i};

        alignment_p = alignments(strcmp(alignments.utter_info, fname_p),:);
        alignment_l = alignments(strcmp(alignments.utter_info, fname_l),:);
        
        audio_p = audioread(fullfile("lombardgrid_paired","audio", [fname_p, '.wav']));
        audio_l = audioread(fullfile("lombardgrid_paired","audio", [fname_l, '.wav']));

        %%Check alignments are not empty and have the same length
        isNotEmpty = (height(alignment_p) ~= 0) && (height(alignment_l) ~=0);
    
        if isNotEmpty
            
            % Obtain index of each segment 
            utterStart_p = alignment_p.offset(1);
            utterEnd_p = alignment_p.offset(end) + alignment_p.duration(end) - 1;
            utterStart_l = alignment_l.offset(1);
            utterEnd_l = alignment_l.offset(1) + alignment_l.duration(end) - 1;

            
            if (utterEnd_p <= length(audio_p)) && (utterEnd_l <= length(audio_l))...
               && ((utterStart_p > 0) && (utterStart_l > 0))

                % Add 30ms silence to the edges of the audio signal to do windowing 
                audio_p = [zeros(edgeFrame, 1); audio_p; zeros(edgeFrame, 1)];
                audio_l = [zeros(edgeFrame, 1); audio_l; zeros(edgeFrame, 1)];
                audio_m = audio_p;

                utterStart_p = utterStart_p + edgeFrame;
                utterEnd_p =  utterEnd_p + edgeFrame;
                utterStart_l = utterStart_l + edgeFrame;
                utterEnd_l =  utterEnd_l + edgeFrame;
                    
                    
                utter_p = audio_p( (utterStart_p - edgeFrame * 0.5) : (utterEnd_p + edgeFrame * 0.5) );  
                utter_l = audio_l( (utterStart_l - edgeFrame * 0.5) : (utterEnd_l + edgeFrame * 0.5) );

                utter_p = edgeWindowing(utter_p, edgeFrame * 2);
                utter_l = edgeWindowing(utter_l, edgeFrame * 2);

                % Measure the spectral-tilts of plain/lombard segment
                tilt_p = myGetSpectralTilt(utter_p, fs, method);
                tilt_l = myGetSpectralTilt(utter_l, fs, method);

                if strcmpi(lpc,'lpc')
                    
                    utter_m = LPCTransplantation(utter_p, utter_l, emp_coeff, fs);
                else

                    utter_m = iterateSTfilter(utter_p, fs, tilt_p, tilt_l, method);
                end


                tilt_m = myGetSpectralTilt(utter_m, fs, method);

                % Update audio_m

                audio_m((utterStart_p - edgeFrame * 0.5) : (utterEnd_p + edgeFrame * 0.5)) = ...
                edgeWindowing(audio_m((utterStart_p - edgeFrame * 0.5) : (utterEnd_p + edgeFrame * 0.5)), edgeFrame, 'inverse');

                audio_m((utterStart_p - edgeFrame * 0.5) : (utterEnd_p + edgeFrame * 0.5)) =...
                audio_m((utterStart_p - edgeFrame * 0.5) : (utterEnd_p + edgeFrame * 0.5)) + utter_m;
                      
                
                % Delete the silences added on the edges of audio
                audio_m = audio_m( edgeFrame : end - edgeFrame );

                % Gain match
                audio_m = audio_m * rms(audio_p)/rms(audio_m);
        
                
                tilt_err_utter = tilt_l - tilt_m;

                result.utterance(i) = fname_p;
                result.tilt_p(i) = tilt_p;
                result.tilt_l(i) = tilt_l;
                result.tilt_m(i) = tilt_m;
                result.abserr_utter(i) = tilt_err_utter;
    
                % Store modified utterance
                audioName = fullfile(savefolder, [speaker '_m_' utterance '.wav']);
                audiowrite(audioName, audio_m, fs);
    
            else
                error(errorIndex,:) = [fname_p -1];
                errorIndex = errorIndex + 1;
                disp('error(-1): The index of the last segment (segEnd(end)) exceeds the length of the audio');
            end
    
        else
            error(errorIndex,:) = [fname_p -2];
            errorIndex = errorIndex + 1;
            disp('error(-2): the alignments do not have the same length or do not have any data');
        end
    
    end
    
    % remove empty rows
    result = rmmissing(result);
    error = rmmissing(error);
    save(fullfile(savefolder,[savefile, '.mat']),'result','error','-mat');

end
