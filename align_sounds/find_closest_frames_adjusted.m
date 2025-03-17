function [frames,difference2] = find_closest_frames_adjusted(frame_times,data,min_distance)
frames = zeros(1,length(data));
count = 0;
for datapoint = data
    count = count+1;
    difference = frame_times - datapoint;
    temp = find(difference>0,1);
    difference2 = difference(temp);
    if frame_times(temp) - datapoint < min_distance
        frames(count) = temp;
    else
        frames(count) =nan;
    end
end