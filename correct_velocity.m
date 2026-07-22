function [corrected_velocity]=correct_velocity(raw_velocity,calibration_info)
% function that corrects for offset and gain of the velocity

% INPUTS:
% raw_velocity: 3 x total frames array of uncorrected velocity that comes straight from the digidata or wavesurfer file 
% calibration_info: structure containing the gain and offset correction for a particular microscope for a particular time period

% OUTPUTS: 
% final_velocity: ([roll pitch yaw], time) array to be used in analysis 

%% Make Variable
corrected_velocity=nan(size(raw_velocity,1),size(raw_velocity,2)); 


%% Correct Offset

raw_velocity(1,:)=raw_velocity(1,:) + calibration_info.pitchoffset; 
raw_velocity(2,:)=raw_velocity(2,:) + calibration_info.rolloffset;
raw_velocity(3,:)=raw_velocity(3,:) + calibration_info.yawoffset; 

%% Correct Gain 

corrected_velocity(1,:)=raw_velocity(1,:) * calibration_info.pitchgain; 
corrected_velocity(2,:)=raw_velocity(2,:) * calibration_info.rollgain;
corrected_velocity(3,:)=raw_velocity(3,:) * calibration_info.yawgain; 


end
