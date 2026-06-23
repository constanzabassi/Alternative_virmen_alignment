function [all_celltypes,allcells_st,dff_st,deconv_st,deconv_st_nogap,stim_info,allcells_st_active,dff_st_active,deconv_st2,deconv_st_nogap2,stim_info2,active,pass,dff_st_noise,deconv_st_noise,deconv_st_nogap_noise,stim_info_noise] = pool_activity_sounds_simple(mouse_date,server,before_after_frames)
dff_st = {};
allcells_st=[];
dff_st_active = {};
allcells_st_active=[];



for m = 1:length(mouse_date)
    mm = mouse_date(m)
    mm = mm{1,1};
    ss = server(m);
    ss = ss {1,1};
%     load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/bad_frames.mat'));
    deconv = load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/deconv/deconv.mat'));
    load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/dff.mat'));
    deconv = deconv.deconv;

    %load trial stuff
    load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/VR/vr_sound_frames.mat'));
    load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/VR/imaging.mat')); %to get trials that actually have VR
    passive_frames = load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/passive/passive_frames.mat')).passive_frames;

    %load alignment info
    load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/alignment_info.mat'));
    passive_temp_dir = cellfun(@(x) contains(x,'passive'),{alignment_info.sync_id},'UniformOutput',false);
    passive_dir = find([passive_temp_dir{1,:}]);

    vr_sound_frames_updated = fix_vr_sound_frames(vr_sound_frames, imaging, alignment_info)

    frame_lengths = [];
    frame_lengths = cellfun(@length,{alignment_info.frame_times}); %across all imaged files 
    frame_lengths = [0,cumsum(frame_lengths)];
    passive_to_add = frame_lengths(passive_dir(1));

    % pass.resp = resp_tr;
    pass.alignment_frames_all = passive_frames.corr_frames;%+passive_to_add;
    pass.trials =  passive_frames.trial_num;
    pass.condition = passive_frames.condition;
    [~,pass.trials_first_repeat] =  unique(passive_frames.trial_num);
    pass.loc_trial = find_trial_conditions(passive_frames.condition,'first_repeat',pass.trials_first_repeat);
    

    % active.resp = resp_tr;
    
    active.alignment_frames_all = vr_sound_frames_updated.corr_frames;
    active.trials =  vr_sound_frames_updated.trial_num; %using trials that have VR data
    [~,active.trials_first_repeat] =  unique(active.trials);
    active.loc_trial = find_trial_conditions(vr_sound_frames_updated.condition,'first_repeat',active.trials_first_repeat);
    

    %load passive variable
% 
% pass.resp = resp_tr;
% pass.loc_trial = loc_trial;
% pass.ctrl_output = control_output_all;
% pass.sound_onsets_all = sound_onsets_all;
% pass.alignment_frames_all = alignment_frames_all;
% 
% %load active variables
% 
% active.resp = resp_tr;
% active.loc_trial = loc_trial;
% active.ctrl_output = control_output_all;
% active.sound_onsets_all = sound_onsets_all;
% active.alignment_frames_all = alignment_frames_all;


    %do dff first!
    %passive stim is sound 1, ctrl is sound 2! for pool activity sound
    %here I am using stim as sound+noise trial and ctrl as sound alone
%     [allcells,allcells_nogap,alignment_frames] = optoalign_function(pass.sound_onsets_all{1,m}(pass.loc_trial{m,1})', pass.sound_onsets_all{1,m}(pass.loc_trial{m,2})', pass.alignment_frames_all,dff,deconv, before_after_frames(1),before_after_frames(2));
    [allcells,allcells_nogap,alignment_frames] = optoalign_function([pass.loc_trial{3}',pass.loc_trial{4}'], [pass.loc_trial{1}',pass.loc_trial{2}'], [pass.alignment_frames_all(:,1),pass.alignment_frames_all(:,1)],dff,deconv, before_after_frames(1),before_after_frames(2),[0,0]);
    [stim_matrix, ctrl_matrix,z_stim_matrix, z_ctrl_matrix] = make_tr_cel_time(allcells,1); %1 is dff

    dff_st{m}.stim = stim_matrix;
    dff_st{m}.ctrl = ctrl_matrix; 
    dff_st{m}.z_stim = z_stim_matrix;
    dff_st{m}.z_ctrl =  z_ctrl_matrix;
    stim_info{m,1} = pass.alignment_frames_all; %bad_frames
    stim_info{m,2} = [pass.loc_trial{3}',pass.loc_trial{4}']; %exp
    stim_info{m,3} = [pass.loc_trial{1}',pass.loc_trial{2}']; %nonexp
    allcells_st= [allcells_st,allcells];
    [deconv_stim, deconv_control] = make_tr_cel_time(allcells_nogap,0);

    deconv_st_nogap{m}.stim = deconv_stim;
    deconv_st_nogap{m}.ctrl = deconv_control;

    %do the same but now interpolate?
    [deconv_stim2, deconv_control2] = make_tr_cel_time(allcells,0);
    deconv_st{m}.stim = deconv_stim2;
    deconv_st{m}.ctrl = deconv_control2;

    %active! stim is sound 1, ctrl is sound 2!
    [allcells,allcells_nogap,alignment_frames] = optoalign_function([active.loc_trial{3}',active.loc_trial{4}'], [active.loc_trial{1}',active.loc_trial{2}'], [active.alignment_frames_all(:,1),active.alignment_frames_all(:,1)],dff,deconv, before_after_frames(1),before_after_frames(2),[0,0]);
    [stim_matrix, ctrl_matrix,z_stim_matrix, z_ctrl_matrix] = make_tr_cel_time(allcells,1); %1 is dff

    dff_st_active{m}.stim = stim_matrix;
    dff_st_active{m}.ctrl = ctrl_matrix; 
    dff_st_active{m}.z_stim = z_stim_matrix;
    dff_st_active{m}.z_ctrl =  z_ctrl_matrix;
    stim_info2{m,1} = active.alignment_frames_all; %bad_frames
    stim_info2{m,2} = [active.loc_trial{3}',active.loc_trial{4}']; %exp
    stim_info2{m,3} = [active.loc_trial{1}',active.loc_trial{2}']; %nonexp
    allcells_st_active= [allcells_st_active,allcells];

    [deconv_stim, deconv_control] = make_tr_cel_time(allcells_nogap,0);

    deconv_st_nogap2{m}.stim = deconv_stim;
    deconv_st_nogap2{m}.ctrl = deconv_control;

    %do the same but now interpolate?
    [deconv_stim2, deconv_control2] = make_tr_cel_time(allcells,0);
    deconv_st2{m}.stim = deconv_stim2;
    deconv_st2{m}.ctrl = deconv_control2;

    all_celltypes{m}.pyr_cells = 1:size(dff,1);


    %get noise only trials
    spont_alignment = pass.alignment_frames_all(:,1);
    spont_alignment(pass.loc_trial{1}) = pass.alignment_frames_all(pass.loc_trial{1},1)+100; %align to nothing (no sound)
    [allcells_noise,allcells_nogap_noise,alignment_frames_noise] = optoalign_function([pass.loc_trial{5}'], [pass.loc_trial{1}'], [spont_alignment,spont_alignment],dff,deconv, before_after_frames(1),before_after_frames(2),[0,0]);
    [stim_matrix_noise, ctrl_matrix_noise,z_stim_matrix_noise, z_ctrl_matrix_noise] = make_tr_cel_time(allcells_noise,1); %1 is dff


    dff_st_noise{m}.stim = stim_matrix_noise;
    dff_st_noise{m}.ctrl = ctrl_matrix_noise; 
    dff_st_noise{m}.z_stim = z_stim_matrix_noise;
    dff_st_noise{m}.z_ctrl =  z_ctrl_matrix_noise;
    stim_info_noise{m,1} = pass.alignment_frames_all; %bad_frames
    stim_info_noise{m,2} = [pass.loc_trial{5}']; %exp
    stim_info_noise{m,3} = [pass.loc_trial{1}']; %nonexp
    allcells_st= [allcells_st,allcells_noise];
    [deconv_stim, deconv_control] = make_tr_cel_time(allcells_nogap_noise,0);

    deconv_st_nogap_noise{m}.stim = deconv_stim;
    deconv_st_nogap_noise{m}.ctrl = deconv_control;

    %do the same but now interpolate?
    [deconv_stim2, deconv_control2] = make_tr_cel_time(allcells_noise,0);
    deconv_st_noise{m}.stim = deconv_stim2;
    deconv_st_noise{m}.ctrl = deconv_control2;

end
%load red cells
% all_celltypes = {};
% for m = 1:length(mouse_date)
%     mm = mouse_date(m);
%     mm = mm{1,1};
%     ss = server(m);
%     ss = ss {1,1};
%     if isdir(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/red_variables/'))==1
%         load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/red_variables/pyr_cells.mat'));
%         load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/red_variables/tdtom_cells.mat'));
%         load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/red_variables/mcherry_cells.mat'));
%         
%         all_celltypes{m}.pyr_cells = pyr_cells';
%         all_celltypes{m}.som_cells= mcherry_cells;
%         all_celltypes{m}.pv_cells= tdtom_cells;
%     
%         total_sum = [length(all_celltypes{m}.som_cells)+length(all_celltypes{m}.pyr_cells)+length(all_celltypes{m}.pv_cells)];
%         if total_sum == size(deconv_st{m}.stim,2) && total_sum == size(dff_st{m}.stim,2)
%             fprintf([num2str(mm) ': cell numbers are a match!\n'])
%         else
%             fprintf([num2str(mm) ': cell numbers dont match!\n'])
%         end
%     end

end