function [warning_index,catch_index, go_index,startle_index,warning_go, startle_sig] = sound_markers_siren(sound_sig1,sound_sig2,fs,cue)
%This code outputs the warning, go and startle cues
% Inputs: The sound signal and the sampling frequency
% Outputs: the locations for warning, go and startle cues

%% Loading the data

% [filename,pathname] = uigetfile ('*.bdf','Please select the file to be analyzed.');
% fileSlash = GetFileSlash(pathname);
% headerFile=strcat(pathname,fileSlash,filename);
% % 
% % nEMGs = 12; %change depending on the number of muscles recorded
% % [data,numChan,labels,txt,fs,gain,prefiltering,ChanDim] = eeg_read_bdf(filename,'all','n');
% % EMG_data = data(1:nEMGs*2, :); %*2 because there are two electrodes for one muscle
% sound_index1 =  find(strcmp(labels, 'Ana6')); %analog 6 is what we use to record the sound
% sound_sig1= data(sound_index1,:);% sound signal 
% sound_index2 =  find(strcmp(labels, 'Ana7')); %analog 6 is what we use to record the sound
% sound_sig2= data(sound_index2,:);% sound signal 
%% Find startle cues
%Startle
% startle_coef = 0.1;
% startle_baseLine = max(sound_sig2(1:fs));
% startle_sound = sound_sig2 - startle_baseLine;
% plot(startle_sound)
% startle_threshold = abs(mean(startle_sound));%startle_coef*mean(startle_sound);
% %[~, index1] = findpeaks(startle_sound, 'MinPeakHeight',startle_threshold);
% [~, index1] = find(startle_sound > startle_threshold);
% %startle_index = find(sound_sig > startle_threshold);
% %index2 = startle_sound(index1); %> startle_threshold;
% startle_index = index1;%(index2);
% indices = find((diff(startle_index) > 0.5*fs) ==1);
% startle_index = [startle_index(1) startle_index(indices+1)];

%trigger
startle_baseLine = max(sound_sig2(fs*3:4*fs));
startle_sound = sound_sig2(fs*3:end)-startle_baseLine;
plot(startle_sound)
startle_threshold = abs(mean(startle_sound(fs:2*fs)));%startle_coef*mean(startle_sound);
%[~, index1] = findpeaks(startle_sound, 'MinPeakHeight',startle_threshold);
[~, index1] = find(abs(startle_sound) > startle_threshold);
%startle_index = find(sound_sig > startle_threshold);
%index2 = startle_sound(index1); %> startle_threshold;
startle_index = index1+((fs*3)-1);%(index2);
indices = find((diff(startle_index) > 0.5*fs) ==1);
startle_index = [startle_index(1) startle_index(indices+1)];
plot(startle_sound)
hold on
plot(startle_index-(fs*3), startle_sound(startle_index-(fs*3)), '*')

%% Find warning and go cues
%Removing all startle cues
go_baseLine = mean(sound_sig1(1:1*fs));
warning_go = sound_sig1;
cue_duration = 0.065*fs; %sound cue lasts for 40ms but because of noise/delay let's use 0.06
% startle_cue_duration = 0.7*fs; %startle sound cue lasts for 60ms 
% for i=1:length(startle_index)
%     warning_go(fix(startle_index(i)-(startle_cue_duration/2)):fix(startle_index(i)+((3/2) *startle_cue_duration)))= ...
%     sound_sig1(fix(startle_index(i)-(4*startle_cue_duration)):fix(startle_index(i)-(2*startle_cue_duration)));
% end
figure(2)
warning_go = warning_go-go_baseLine;
plot(warning_go)

%mean_sig = mean(warning_go(1:round(2*cue_duration)));
%std_sig = std(warning_go(1:round(2*cue_duration)));
%wg_coef = 1;
%warning_go_thr = mean_sig + (wg_coef*std_sig);

lowerbound = 3*fs; %shortest time for a second sound
%find peaks of the signal
peakProm = 1800; %change if baseline of noise, is more or less
[~, wg_index1] = findpeaks(warning_go,"MinPeakDistance",lowerbound, ...
                "MinPeakProminence",peakProm);
%wg_index2 = warning_go(wg_index1) > warning_go_thr;
%wg_index = wg_index1(wg_index2);
%wg_index_diff1 = diff(wg_index);

%lowerbound = 1*fs; %shortest time for a second sound
%index_fault = find(wg_index_diff1<lowerbound)+1;
%wg_index(index_fault) = []; %delete those that are too close
wg_index = wg_index1;
wg_index_diff2 = diff(wg_index);
index_main = sort([wg_index, startle_index]);
index_diff = diff(index_main);

upperbound = 6*fs; %adding .5 because of delay interval is 3.5 - 4.5
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
%amplitude = 1/100;
startle_sig = sound_sig2-startle_baseLine;%*amplitude;
% sound_signal = warning_go+scaled_startle_sig;
sas_baselineCorrect = (sound_sig2 - startle_baseLine) + 10000;

% for i=1:length(startle_index)
%     scaled_sound_sig(fix(startle_index(i)-(cue_duration/2)):fix(startle_index(i)+((3/2) *cue_duration)))=(sound_sig1(fix(startle_index(i)-(cue_duration/2)):fix(startle_index(i)+((3/2) *cue_duration))))*amplitude;
% end
% plot(startle_sound)
% hold on
% plot(scaled_sound_sig)
if cue ~= 5
    figure(123)
    plot(warning_go,'k')
    hold on
    plot(sas_baselineCorrect,'g')
    hold on
    a(1) = plot(go_index, warning_go(go_index),'*g', 'DisplayName', 'Go');
    a(2) = plot(warning_index, warning_go(warning_index),'*k', 'DisplayName', 'Warning');
    a(3) = plot(startle_index, sas_baselineCorrect(startle_index), '*r','DisplayName', 'Startle');
    legend(a)
    
end


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