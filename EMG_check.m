function EMG_check(filename, headerFile, filt_outputData,~,filt_outputData_rect,catch_marker_index,go_marker_index,startle_marker_index ,~, fs, cue, EMGsOrder) 
clc
% if ~isempty(emg_mvc)
%     mvcs = ([emg_mvc; 0; 0])';
% else
%     disp('Warning: You are processing this data without the MVCs')
%     mvcs = '';
% end

%% Finding the EMG onsets and offsets  

maxCatchEMGonset_duration = 3.5*fs; %this is the 3000 - 8000ms *fs
maxEMGonset_duration = 2*fs; %this is the 5000 - 8000ms *fs
maxEMG_duration = 20*fs; %this is the 100ms *fs 
emg_preLength = 1*fs; %1 %500ms*fs 
scm_pre_length = 0.05*fs; %20ms*fs


% Self Initiated
maxRT = 150; %max RT = 150 for startle 
if cue == 1 
    sound_cue = 'SI';
    manual_onsetCheck_CCC(headerFile,filename,filt_outputData,filt_outputData_rect,catch_marker_index, ...
        maxCatchEMGonset_duration,maxEMG_duration,EMGsOrder,sound_cue,0,scm_pre_length,fs);
    
%Go trials
elseif cue == 2 
    sound_cue = 'Go';
    [scm_EMGs,scm_threshold,SCM_RT,RT] = manual_onsetCheck(headerFile,filename,filt_outputData,filt_outputData_rect,go_marker_index, ...
        emg_preLength,maxEMGonset_duration,EMGsOrder,sound_cue,emg_preLength,scm_pre_length,fs);
    startleResponse_check(filename,headerFile, scm_EMGs, SCM_RT, scm_threshold, RT, maxRT,fs,sound_cue);

%Startle Trials
elseif cue == 3 
    sound_cue = 'Startle';
    [scm_EMGs,scm_threshold,SCM_RT,RT] = manual_onsetCheck(headerFile,filename,filt_outputData,filt_outputData_rect,startle_marker_index, ...
        emg_preLength,maxEMGonset_duration,EMGsOrder,sound_cue,emg_preLength,scm_pre_length,fs);
    startleResponse_check(filename,headerFile, scm_EMGs, SCM_RT, scm_threshold, RT, maxRT,fs,sound_cue);
    
%Go&Startle Trials
elseif cue == 4 
    sound_cue = 'Go';
    [scm_EMGs,scm_threshold,SCM_RT,RT] = manual_onsetCheck(headerFile,filename,filt_outputData,filt_outputData_rect,go_marker_index, ...
        emg_preLength,maxEMGonset_duration,EMGsOrder,sound_cue,emg_preLength,scm_pre_length,fs);
    startleResponse_check(filename,headerFile, scm_EMGs, SCM_RT, scm_threshold, RT, maxRT,fs,sound_cue);

    sound_cue = 'Startle';
    [scm_EMGs,scm_threshold,SCM_RT,RT] = manual_onsetCheck(headerFile,filename,filt_outputData,filt_outputData_rect,startle_marker_index, ...
        emg_preLength,maxEMGonset_duration,EMGsOrder,sound_cue,emg_preLength,scm_pre_length,fs);
    startleResponse_check(filename,headerFile, scm_EMGs, SCM_RT, scm_threshold, RT, maxRT,fs,sound_cue);

elseif cue == 5 
    sound_cue = 'ClassicSAS';
    classicSAS_control(headerFile,filename,filt_outputData,filt_outputData_rect,startle_marker_index, ...
        emg_preLength,maxEMGonset_duration,EMGsOrder,sound_cue,scm_pre_length,fs);
    %startleResponse_check(filename,headerFile, scm_EMGs, SCM_RT, scm_threshold, RT, maxRT,fs,sound_cue);
      
end


