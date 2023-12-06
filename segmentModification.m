function segmentModification(corpus, alignments, savefolder, savefile, method, lpc)
% input
% corpus: using corpus
% maxRelErr: the target modifiction relative error
% savefolder: 
% savefile: 

    if nargin < 6, lpc = ''; end

    % coefficient for pre/de emphasis filter
    emp_coeff = -0.97; 
    fs = 16000;
    
    result = table('Size', [height(corpus) * 30 11], ...
                    'VariableTypes', {'string','string','double','double','double','double','double','double','double','double','double'}, ...
                    'VariableNames',{'utterance','phone','tilt_p','tilt_l','tilt_m','abserr_seg','abserr_utter','segStart_p','segEnd_p','segStart_l','segEnd_l'});
    
    error = table('Size', [height(corpus) 2], 'VariableTypes',{'string','double'}, ...
                                               'VariableNames',{'utterance','error_type'});

    size = height(corpus);
    
    mkdir(savefolder);
    
    utter_start = 1;
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

        if height(alignment_p) == 0 || height(alignment_l) == 0
            
            plainLength = 0;
            lombardLength = 0;
    
        else

            plainLength = height(alignment_p);
            lombardLength = height(alignment_l);
    
        end
    
        %%Check alignments are not empty and have the same length
        isNotEmpty = (plainLength ~= 0) && (lombardLength ~=0);
        isSameLength = (plainLength == lombardLength); 
    
        if isNotEmpty  && isSameLength
            
            % Obtain index of each segment 
            segStart_p = alignment_p.offset;
            segEnd_p = alignment_p.offset + alignment_p.duration - 1;
            segStart_l = alignment_l.offset;
            segEnd_l = alignment_l.offset + alignment_l.duration - 1;

            
            if (segEnd_p(end) <= length(audio_p)) && (segEnd_l(end) <= length(audio_l))...
               && ((segStart_p(1) > 0) && (segStart_l(1) > 0))

                tilt_p = cell(plainLength,1);
                tilt_l = tilt_p;
                tilt_m = tilt_p;

                % Add 30ms silence to the edges of the audio signal to do windowing 
                audio_p = [zeros(edgeFrame, 1); audio_p; zeros(edgeFrame, 1)];
                audio_l = [zeros(edgeFrame, 1); audio_l; zeros(edgeFrame, 1)];
                audio_m = audio_p;

                segStart_p = segStart_p + edgeFrame;
                segEnd_p =  segEnd_p + edgeFrame;
                segStart_l = segStart_l + edgeFrame;
                segEnd_l =  segEnd_l + edgeFrame;
                    
                for j = 1 : plainLength
                    
                    seg_p = audio_p( (segStart_p(j) - edgeFrame * 0.5) : (segEnd_p(j) + edgeFrame * 0.5) );  
                    seg_l = audio_l( (segStart_l(j) - edgeFrame * 0.5) : (segEnd_l(j) + edgeFrame * 0.5) );

                    seg_p = edgeWindowing(seg_p, edgeFrame * 2);
                    seg_l = edgeWindowing(seg_l, edgeFrame * 2);

                    % Measure the spectral-tilts of plain/lombard segment
                    tilt_p{j} = myGetSpectralTilt(seg_p, fs, method);
                    tilt_l{j} = myGetSpectralTilt(seg_l, fs, method);


                    if strcmpi(lpc,'lpc')
                        
                        seg_m = LPCTransplantation(seg_p, seg_l, emp_coeff, fs);
                    else

                        seg_m = iterateSTfilter(seg_p, fs, tilt_p{j}, tilt_l{j}, method);
                    end

                    tilt_m{j} = myGetSpectralTilt(seg_m, fs, method);

                    % Update audio_m

                    audio_m((segStart_p(j) - edgeFrame * 0.5) : (segEnd_p(j) + edgeFrame * 0.5)) = ...
                    edgeWindowing(audio_m((segStart_p(j) - edgeFrame * 0.5) : (segEnd_p(j) + edgeFrame * 0.5)), edgeFrame, 'inverse');

                    audio_m((segStart_p(j) - edgeFrame * 0.5) : (segEnd_p(j) + edgeFrame * 0.5)) =...
                    audio_m((segStart_p(j) - edgeFrame * 0.5) : (segEnd_p(j) + edgeFrame * 0.5)) + seg_m;
                      
                end
                
                % Delete the silences added on the edges of audio
                audio_m = audio_m( edgeFrame : end - edgeFrame );

                % Gain match
                audio_m = audio_m * rms(audio_p)/rms(audio_m);
        
                utter_end = utter_start + plainLength - 1;
        
                result(utter_start : utter_end,{'utterance','phone'}) = alignment_p(:, {'utter_info','phone'});
                result(utter_start : utter_end,{'tilt_p','tilt_l','tilt_m'}) = [tilt_p tilt_l tilt_m];
                result(utter_start : utter_end,{'segStart_p','segEnd_p'}) = num2cell([segStart_p segEnd_p]);
                result(utter_start : utter_end,{'segStart_l','segEnd_l'}) = num2cell([segStart_l segEnd_l]);
    
    
                tilt_err_seg = cell2mat(tilt_l) - cell2mat(tilt_m);
                result{utter_start:utter_end, 'abserr_seg'} = tilt_err_seg;
                
                result{utter_start:utter_end, 'abserr_utter'} = mean(abs(tilt_err_seg));
            
                % Store modified utterance
                audioName = fullfile(savefolder, [speaker '_m_' utterance '.wav']);
                audiowrite(audioName, audio_m, fs);
        
                utter_start = utter_end + 1;
    
            else
                error(errorIndex,:) = [corpus.FNAME_P(i) -1];
                errorIndex = errorIndex + 1;
                disp('error(-1): The index of the last segment (segEnd(end)) exceeds the length of the audio');

                disp('plain');
                disp(segStart_p(1))
                disp(length(audio_p) - segEnd_p(end))
                disp('lombard');
                disp(segStart_l(1))
                disp(length(audio_l) - segEnd_l(end))
            end
    
        else
            error(errorIndex,:) = [corpus.FNAME_P(i) -2];
            errorIndex = errorIndex + 1;
            disp('error(-2): the alignments do not have the same length or do not have any data');
        end
    
    end
    
    % remove empty rows
    result = rmmissing(result);
    error = rmmissing(error);
    save(fullfile(savefolder,[savefile, '.mat']),'result','error','-mat');

end
