
function [alignment_p, alignment_l, audio_p, audio_l, Fs] = get_audio_and_alignment(speaker, utterance)
% This function returns alignments and plain/lombard speeches that are produced by the same speaker.
% Alignment is the table that contains the duration of the segment, the offset, and the phone
% 
% INPUT  speaker:
%       utterance:

wildcard = ['**/' speaker '*' utterance '*'];

files_info = dir(wildcard);

if length(files_info) == 0
    disp([speaker '_'  utterance ' does not exist.' ]);
    return
end

alignment_p = table([],[],[],[],'VariableNames', ...
            {'utter_info','duration','offset','phone'});
alignment_l = alignment_p;

audio_p = [];
audio_l = [];
Fs = 0;


voiceless = {'SIL','SIL_S'};
% voiceless = {'ch','f','k','p','s','sh','t','th','SIL','SIL_S'};
IsVoiced = @(alignment) ~contains(alignment.phone, voiceless);


for n = 1:length(files_info)

    fullpath = [files_info(n).folder '\' files_info(n).name];

    if contains(fullpath, 'json') 
        
        fileID = fopen(fullpath);
        raw = fread(fileID, inf);
        str = char(raw');
        fclose(fileID);

        if contains(fullpath, '_p_')
            alignment_p = struct2table(jsondecode(str).( [speaker '_p_' utterance] ));

            %%delete voiceless and silence from the alignments
            alignment_p = alignment_p(IsVoiced(alignment_p), :);

            utter_info = table(repmat({[speaker '_p_' utterance]}, height(alignment_p), 1) ...
                    ,'VariableNames', {'utter_info'});

            alignment_p = [utter_info alignment_p(:,{'phone','duration','offset'})];
        else
            alignment_l = struct2table(jsondecode(str).( [speaker '_l_' utterance] ));
            
            %%delete voiceless and silence from the alignments
            alignment_l = alignment_l(IsVoiced(alignment_l), :);
            
            utter_info = table(repmat({[speaker '_l_' utterance]}, height(alignment_l), 1) ...
                    ,'VariableNames', {'utter_info'});

            alignment_l = [utter_info alignment_l(:,{'phone','duration','offset'})];
        end

    else
        
        if contains(fullpath, '_p_')
            [audio_p, Fs] = audioread(fullpath);
        else
            [audio_l, Fs] = audioread(fullpath);
        end

    end
     
end

end