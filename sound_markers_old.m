function [warning_index,catch_index, go_index,startle_index, sound_sig] = sound_markers_old(sound_sig,fs)
%This code outputs the warning, go and startle cues
% Inputs: The sound signal and the sampling frequency
% Outputs: the locations for warning, go and startle cues

%% Loading the data

% [filename,pathname] = uigetfile ('*.bdf','Please select the file to be analyzed.');
% fileSlash = GetFileSlash(pathname);
% headerFile=strcat(pathname,fileSlash,filename);
% 
% nEMGs = 12; %change depending on the number of muscles recorded
% [data,numChan,labels,txt,fs,gain,prefiltering,ChanDim] = eeg_read_bdf(filename,'all','n');
% EMG_data = data(1:nEMGs*2, :); %*2 because there are two electrodes for one muscle
% sound_index =  find(strcmp(labels, 'Ana6')); %analog 6 is what we use to record the sound
% sound_sig= data(sound_index,:);% sound signal 

%% Find startle cues
%Startle
startle_coef = 0.8;
startle_threshold = abs(startle_coef*max(abs(sound_sig)));
%[~, index1] = findpeaks(sound_sig);
indices = find(abs(sound_sig) > startle_threshold);
index= find((diff(indices) > 0.5*fs) ==1);
startle_index = [indices(1) indices(index+1)];
plot(sound_sig)
hold
plot(startle_index, sound_sig(startle_index), '*')

%startle_index = find(sound_sig > startle_threshold);
% index2 = sound_sig(index1) > startle_threshold;
% startle_index = index1(index2);
% 
% startle_baseLine = max(sound_sig2(fs*3:4*fs));
% startle_sound = sound_sig2(fs*3:end)-startle_baseLine;
% plot(startle_sound)
% startle_threshold = 100*abs(mean(startle_sound(fs:2*fs)));%startle_coef*mean(startle_sound);
% %[~, index1] = findpeaks(startle_sound, 'MinPeakHeight',startle_threshold);
% [~, index1] = find(abs(startle_sound) > startle_threshold);
% %startle_index = find(sound_sig > startle_threshold);
% %index2 = startle_sound(index1); %> startle_threshold;
% startle_index = index1+((fs*3)-1);%(index2);
% indices = find((diff(startle_index) > 0.5*fs) ==1);
% startle_index = [startle_index(1) startle_index(indices+1)];
% plot(startle_sound)
% hold on
% plot(startle_index-(fs*3), startle_sound(startle_index-(fs*3)), '*')

%% Find warning and go cues
%Removing all startle cues
warning_go = sound_sig;
cue_duration = 0.06*fs; %sound cue lasts for 40ms but because of noise/delay let's use 0.06
for i=1:length(startle_index)
    warning_go(fix(startle_index(i)-(cue_duration/2)):fix(startle_index(i)+((3/2) *cue_duration)))=sound_sig(fix(startle_index(i)-(4*cue_duration)):fix(startle_index(i)-(2*cue_duration)));
end
plot(warning_go)

mean_sig = mean(warning_go(1:round(2*cue_duration)));
std_sig = std(warning_go(1:round(2*cue_duration)));
wg_coef = 4;
warning_go_thr = mean_sig + (wg_coef*std_sig);

[~, wg_index1] = findpeaks(warning_go);
wg_index2 = warning_go(wg_index1) > warning_go_thr;
wg_index = wg_index1(wg_index2);
wg_index_diff1 = diff(wg_index);

lowerbound = 1*fs; %shortest time for a second sound
index_fault = find(wg_index_diff1<lowerbound)+1;
wg_index(index_fault) = []; %delete those that are too close

wg_index_diff2 = diff(wg_index);
index_main = sort([wg_index, startle_index]);
index_diff = diff(index_main);

upperbound = 4.2*fs; %adding .5 because of delay interval is 1-3
warning_index = wg_index;
go = find(discretize(wg_index_diff2, [lowerbound, upperbound])==1)+1;
index_3 = discretize(index_diff, [lowerbound, upperbound])==1; %warning cues for second sound
go_index = wg_index(go);

warning_index(go) = []; %removing the go cues
catch_index = index_main; 
catch_index(index_3) = [];
index_4 = ismember(catch_index, go_index) | ismember(catch_index, startle_index) ==1;
catch_index(index_4) = [];

%% Plotting
amplitude = 0.00001;
scaled_sound_sig = sound_sig;
for i=1:length(startle_index)
    scaled_sound_sig(fix(startle_index(i)-(cue_duration/2)):fix(startle_index(i)+((3/2) *cue_duration)))=(sound_sig(fix(startle_index(i)-(cue_duration/2)):fix(startle_index(i)+((3/2) *cue_duration))))*amplitude;
end

plot(scaled_sound_sig)
hold on
a(1) = plot(go_index, scaled_sound_sig(go_index),'*g', 'DisplayName', 'Go');
a(2) = plot(warning_index, scaled_sound_sig(warning_index),'*k', 'DisplayName', 'Warning');
a(3) = plot(startle_index, scaled_sound_sig(startle_index), '*r','DisplayName', 'Startle');
legend(a)

% %% Saving
% ind = find (filename == '.');
% outputName = strcat(filename (1:ind-1),'.mat');
% outputData = [EMG_data; sound_sig];
% save (outputName,'outputData','go_index','startle_index','warning_index');

%% Display numbers of each cue
catch_trials = length(catch_index);
fprintf('Number of Startle Trials: %d \n',length(startle_index));
fprintf('Number of Go Trials: %d \n',length(go_index));
fprintf('Number of Catch Trials: %d \n',catch_trials);
fprintf('Number of Trials: %d \n',length(warning_index));

end