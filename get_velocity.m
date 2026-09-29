function[raw_velocity]=get_velocity(alignment_info,calibration_info,sync_base_path,varargin)

%function to return raw_velocity (no gain or offset adjustment)

% TO DO
% - make this compatible with wavesurfer file. 
%%added 3/14/23: optional input ws. If 1, this is a wavesurfer sync

% INPUTS: 
% alignment_info: structure of length T-Series that contains the time of each imaging frame in digidata/wavesurfer coordinates
% calibration_info: can be structure or vector of channel order [pitch, roll, yaw]
% sync_base_path: string containing the file location of wavesurfer/ digidata files 


% OUTPUTS: 
% raw_velocity: array where size(velocity,2) is the number of all frames in all TSeries

%   raw_velocity(1,:)=pitch
%   raw_velocity(2,:)=roll
%   raw_velocity(3,:)=yaw

if length(dir([sync_base_path '*.abf']))>0
    sync_dir = dir([sync_base_path '*.abf']);
    num_syncs = length(sync_dir);
    is_pclamp = 1;
else
    sync_dir = dir([sync_base_path '*.h5']);
    num_syncs = length(sync_dir);
    is_pclamp = 0;
end

%% Get calibration info either from vector or structure 
if isstruct(calibration_info)
    pitch_chan=calibration_info.pitchchannelnum;
    roll_chan=calibration_info.rollchannelnum;
    yaw_chan=calibration_info.yawchannelnum;
else
    pitch_chan=calibration_info(1); 
    roll_chan=calibration_info(2); 
    yaw_chan=calibration_info(3); 
end


%% Load velocity info and put into array 

raw_velocity=[];

for acq_number=1:length(alignment_info)
    
    frame_times=alignment_info(acq_number).frame_times; % list of timepoints in digidata/wavesurfer coordinates where TSeries frame occured
    
    if is_pclamp==0
        data = ws.loadDataFile([sync_base_path alignment_info(acq_number).sync_id]);
        fields = fieldnames(data);
        sweep_id = fields{2};
        sync_data = eval(['data.' sweep_id '.analogScans']);
    else
        [sync_data,~,~] = abfload([sync_base_path alignment_info(acq_number).sync_id]); % load digidata
    end
    velocity=nan(3,length(frame_times)); % make empty velocity array for this TSeries
    
    frame_period=mean(diff(frame_times));% get the periodicity of the frames for the TSeries in digidata coordinates
    
    
    
    for idx=1:length(frame_times)
        
        if frame_times(idx)-round(frame_period/2) < 0
            window=frame_times(idx): frame_times(idx)+round(frame_period/2); %get the window of 1/2 the frame_period on either side of the TSeries frame
        else
        window=frame_times(idx)-round(frame_period/2): frame_times(idx)+round(frame_period/2); %get the window of 1/2 the frame_period on either side of the TSeries frame
        end
        temp_pitch=sync_data(window,pitch_chan);
        temp_roll=sync_data(window,roll_chan);
        temp_yaw=sync_data(window,yaw_chan); 
      
        velocity(1,idx)=mean(temp_pitch);
        velocity(2,idx)=mean(temp_roll);
        velocity(3,idx)=mean(temp_yaw); 
        
    end
    
            
        raw_velocity=cat(2,raw_velocity,velocity); % add velocity vector for this TSeries to larger array 


end


end

