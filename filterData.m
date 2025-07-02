function [filt_outputData, filt_outputData_envelope, filt_outputData_rectify] = filterData(filename, pathname, outputData,catch_marker_index,go_marker_index,startle_marker_index, warning_marker_index, fs, EMGsOrder) 
clc
%if ever you need to run this separately
% [filename, pathname] = uigetfile ('*.mat');
% file2run=strcat(pathname,'/',filename);
% load(file2run)
% fs = 2048;
%% Let's just look at the EMGs too
% hpFilt = designfilt('highpassiir','StopbandFrequency',20, ...
%          'PassbandFrequency',30,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'SampleRate',fs,'DesignMethod','butter');
% % 
% lpFilt = designfilt('lowpassiir','StopbandFrequency',5, ...
%          'PassbandFrequency',3,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'SampleRate',fs,'DesignMethod','butter');

hcut = 20;
lcut = 5;
[bhpFilt ,ahpFilt] = butter(4,hcut*2/fs, 'high');
[blpFilt, alpFilt] = butter(4,lcut*2/fs, 'low');

filt_outputData = zeros(size(outputData));
filt_outputData_envelope = zeros(size(outputData));
filt_outputData_rectify = zeros(size(outputData));

for iChannel = 1:size(outputData,1)
    filt_outputData(iChannel,:) = filtfilt(bhpFilt,ahpFilt, outputData(iChannel,:));

    %% adding this added more noise (low frequency oscillations before emg onset) to my signal
% %     creating notch filters to remove harmonics of 60Hz
%     f0 = [60] %,120,180,240,300]; 
%     for i = f0
%     notchFilt = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',i-1,'HalfPowerFrequency2',i+1, ...
%                'DesignMethod','butter','SampleRate',fs);
%         filt_outputData(iChannel,:) = filtfilt(notchFilt,filt_outputData(iChannel,:));
%     end

    %% lets do the rectification and linear envelope
    filt_outputData_rectify(iChannel,:) = abs(filt_outputData(iChannel,:));
    filt_outputData_envelope(iChannel,:) = filtfilt(blpFilt,alpFilt, filt_outputData_rectify(iChannel,:));
end


%% for debugging purposes
%if you see no emg, just run this to check the emg signal.
%plot_wholeTrial_sasTrigger(filt_outputData,startle_marker_index, go_marker_index,catch_marker_index,EMGsOrder)
%plot_wholeTrial(filt_outputData,startle_marker_index, go_marker_index,catch_marker_index,EMGsOrder)

%% Lets save the filtered data so we don't have to go through all this again
ind = find (filename == '.');
name2save = strcat(pathname, filename (1:ind-1),'_filteredData.mat');
save (name2save,"filt_outputData","filt_outputData_envelope","filt_outputData_rectify", "catch_marker_index","startle_marker_index","go_marker_index","warning_marker_index");



% %%
% N = 2048;
% [bhpFilt ,ahpFilt] = butter(4,hcut*2/fs, 'high');
% [h, w] = freqz(bhpFilt, ahpFilt, N, 'whole',fs);
% 
% % Frequency axis in Hz
% fs = 2048;  % Sampling frequency
% frequencies =  w ;
% magnitude_dB = 20 * log10(abs(h));
% 
% % Find the -3dB point
% [~, idx] = min(abs(magnitude_dB + 3));  % Find index of closest point to -3dB
% frequency_3dB = frequencies(idx);  % Frequency at -3dB
% 
% %
% % Plot frequency response
% figure;
% subplot(2,1,1);
% plot(frequencies, magnitude_dB);
% hold on
% yline(-3, '--r', '-3dB Point');
% xline(frequency_3dB, '--r', 'cut off')
% plot(frequency_3dB, -3, 'ro');  % Mark -3dB point with a red circle
% title('Magnitude Response');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% 
% subplot(2,1,2);
% plot(frequencies, angle(h));
% title('Phase Response');
% xlabel('Frequency (Hz)');
% ylabel('Phase (radians)');
% 
