function updated_condition_array = convert_sound_conditions_simple(condition_array,conditions_per_speaker,speaker_ids)
%speaker_ids = [1,2,4];
total_speakers = 1:4;
updated_condition_array={};
if length(speaker_ids ) == 3
    for t = 1:length(condition_array)
        updated_condition_array{1,t}= condition_array(1,t);
    end
end