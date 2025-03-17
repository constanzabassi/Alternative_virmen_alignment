function [all_celltypes,allcells_st,dff_st,deconv_st,deconv_st_nogap,stim_info,allcells_st2,dff_st2,deconv_st2,deconv_st_nogap2,stim_info2,active,pass] = pool_activity_sounds(mouse_date,server,before_after_frames)
dff_st = {};
allcells_st=[];
dff_st2 = {};
allcells_st2=[];

%load passive variables
load('V:/Connie/results/passive/mod/resp_tr.mat');
load('V:/Connie/results/passive/mod/loc_trial.mat');
load('V:/Connie/results/passive/mod/control_output_all.mat');
load('V:/Connie/results/passive/mod/sound_onsets_all.mat');
load('V:/Connie/results/passive/mod/alignment_frames_all.mat');

pass.resp = resp_tr;
pass.loc_trial = loc_trial;
pass.ctrl_output = control_output_all;
pass.sound_onsets_all = sound_onsets_all;
pass.alignment_frames_all = alignment_frames_all;

%load active variables
load('V:/Connie/results/active/mod/resp_tr.mat');
load('V:/Connie/results/active/mod/loc_trial.mat');
load('V:/Connie/results/active/mod/control_output_all.mat');
load('V:/Connie/results/active/mod/sound_onsets_all.mat');
load('V:/Connie/results/active/mod/alignment_frames_all.mat');

active.resp = resp_tr;
active.loc_trial = loc_trial;
active.ctrl_output = control_output_all;
active.sound_onsets_all = sound_onsets_all;
active.alignment_frames_all = alignment_frames_all;

load('V:\Connie\results\passive\data_info\all_celltypes.mat');

for m = 1:length(mouse_date)
    mm = mouse_date(m)
    mm = mm{1,1};
    ss = server(m);
    ss = ss {1,1};
    load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/bad_frames.mat'));
    deconv = load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/deconv/deconv.mat'));
    load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/dff.mat'));
%     exp = load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/',path_string,'/exp.mat'));
%     nonexp = load(strcat(num2str(ss),'/Connie/ProcessedData/',num2str(mm),'/',path_string,'/nonexp.mat'));
%     field_exp = fields(exp);
%     field_nexp = fields(nonexp);
%     exp = exp.(field_exp{1});
%     nonexp = nonexp.(field_nexp{1});
    deconv = deconv.deconv;

    %do dff first!
    %passive stim is sound 1, ctrl is sound 2!
    [allcells,allcells_nogap,alignment_frames] = optoalign_function(pass.sound_onsets_all{1,m}(pass.loc_trial{m,1})', pass.sound_onsets_all{1,m}(pass.loc_trial{m,2})', pass.alignment_frames_all{1,m},dff,deconv, before_after_frames(1),before_after_frames(2));
    [stim_matrix, ctrl_matrix,z_stim_matrix, z_ctrl_matrix] = make_tr_cel_time(allcells,1); %1 is dff

    dff_st{m}.stim = stim_matrix;
    dff_st{m}.ctrl = ctrl_matrix; 
    dff_st{m}.z_stim = z_stim_matrix;
    dff_st{m}.z_ctrl =  z_ctrl_matrix;
    stim_info{m,1} = pass.alignment_frames_all{1,m}; %bad_frames
    stim_info{m,2} = pass.sound_onsets_all{1,m}(pass.loc_trial{m,1})'; %exp
    stim_info{m,3} = pass.sound_onsets_all{1,m}(pass.loc_trial{m,2})'; %nonexp
    allcells_st= [allcells_st,allcells];
    [deconv_stim, deconv_control] = make_tr_cel_time(allcells_nogap,0);

    deconv_st_nogap{m}.stim = deconv_stim;
    deconv_st_nogap{m}.ctrl = deconv_control;

    %do the same but now interpolate?
    [deconv_stim2, deconv_control2] = make_tr_cel_time(allcells,0);
    deconv_st{m}.stim = deconv_stim2;
    deconv_st{m}.ctrl = deconv_control2;

    %active! stim is sound 1, ctrl is sound 2!
    [allcells,allcells_nogap,alignment_frames] = optoalign_function(active.sound_onsets_all{1,m}(active.loc_trial{m,1})', active.sound_onsets_all{1,m}(active.loc_trial{m,2})', active.alignment_frames_all{1,m},dff,deconv, before_after_frames(1),before_after_frames(2));
    [stim_matrix, ctrl_matrix,z_stim_matrix, z_ctrl_matrix] = make_tr_cel_time(allcells,1); %1 is dff

    dff_st2{m}.stim = stim_matrix;
    dff_st2{m}.ctrl = ctrl_matrix; 
    dff_st2{m}.z_stim = z_stim_matrix;
    dff_st2{m}.z_ctrl =  z_ctrl_matrix;
    stim_info2{m,1} = active.alignment_frames_all{1,m}; %bad_frames
    stim_info2{m,2} = active.sound_onsets_all{1,m}(active.loc_trial{m,1})'; %exp
    stim_info2{m,3} = active.sound_onsets_all{1,m}(active.loc_trial{m,2})'; %nonexp
    allcells_st2= [allcells_st2,allcells];

    [deconv_stim, deconv_control] = make_tr_cel_time(allcells_nogap,0);

    deconv_st_nogap2{m}.stim = deconv_stim;
    deconv_st_nogap2{m}.ctrl = deconv_control;

    %do the same but now interpolate?
    [deconv_stim2, deconv_control2] = make_tr_cel_time(allcells,0);
    deconv_st2{m}.stim = deconv_stim2;
    deconv_st2{m}.ctrl = deconv_control2;

    %organize context_tr!
%     trials_together ={};
%     trials_together{1} = [active.sound_onsets_all{1,m}(active.loc_trial{m,1})',pass.sound_onsets_all{1,m}(pass.loc_trial{m,1})'];
%     trials_together{2} = [active.sound_onsets_all{1,m}(active.loc_trial{m,2})',pass.sound_onsets_all{1,m}(pass.loc_trial{m,2})'];
%     mouse_context_tr{m} = context_tr;


    %save raw data for each mouse
%     mouse{m}.dff = dff;
%     mouse{m}.deconv = deconv;

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