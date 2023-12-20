function plot_random_trials_alignment (imaging,selected_fields)
empty_trials = find(cellfun(@isempty,{imaging.movement_in_imaging_time}));
good_trials =  setdiff(1:length(imaging),empty_trials);
sel_trials = randperm(length(good_trials),9);
selected_trials = good_trials(sel_trials);
figure(1010);clf;
tiledlayout(3,3,"TileSpacing","compact");
for t = 1:9
    imaging_cell = struct2cell(imaging(selected_trials(t)).movement_in_imaging_time);
    nexttile
    for f = 1:length(selected_fields)
        if max(selected_fields(f)) > 1
            plot(rescale(imaging_cell{selected_fields(f),1},0,1));
        else
            plot(imaging_cell{selected_fields(f),1});
        end
    end
end
