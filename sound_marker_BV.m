function [warning_marker_index] = sound_marker_BV(index_markers, sound_sig_main,fs, buffer_gap)
%This code outputs the warning, go and startle cues
% Inputs: The sound signal and the sampling frequency
% Outputs: the locations for warning, go and startle cues
%% first lets filter the signal
[blpFilt, alpFilt] = butter(4,20*2/fs, 'high');
sound_sig = filtfilt(blpFilt,alpFilt, sound_sig_main);
%% Detecting Warning, Go & Startle

soundEnd_time = 10.5*fs; %time from trigger to the end of all sound cues (5.5-10.5)
sound_interval = 0.15*fs; %150ms * fs -> 0.1s *fs (each cue lasts for 40ms)
%(2s for warning, and 1-3.5s for startle/go)
warning_marker = 0;
warning_marker_index = [];
sound_threshold = 0.5e6;
sound_sig_length = length(sound_sig);
n = 1;
for index = index_markers
    buffered_index = index-buffer_gap;
    if buffered_index+soundEnd_time < sound_sig_length
        trigger2cues_interval = [buffered_index:buffered_index+soundEnd_time];
    else
        disp('end of sound signal')
        trigger2cues_interval = [buffered_index:sound_sig_length];
    end
    %detecting the cues
    sound_data_interval = sound_sig(trigger2cues_interval);
    loc = find((sound_sig(trigger2cues_interval))> sound_threshold);
    firstSig_indx = loc(1);
    first_sound = sound_data_interval([firstSig_indx:firstSig_indx+sound_interval]);
 
    %there is no go or startle 
    % Plot to check if the detection is right
    figure 
    subplot(2,1,1)
    plot(sound_data_interval, 'k')
    title(strcat('sound signal in all its glory: ',num2str(n)));
    subplot(2,1,2)
    plot(sound_data_interval, 'k')
    hold on
    plot([firstSig_indx:firstSig_indx+length(first_sound)-1],first_sound, 'r')
    legend('sound signal', 'warning cue')
    title('sound signal with warning cue: ');
   
    warning_marker = warning_marker+1;
    warning_marker_index(1,n) = buffered_index+firstSig_indx;
    
    %selecting onset manually should they not be aligned

    prompt = 'Do you want to reset the onset of the warning cue? 2:yes, 4: no \n';
    answer = input(prompt,'s');
    t = 1;
    while t ==1
        switch answer 
            case '4'
                disp ('You have decided not to reset onset and offset, proceed to the next \n');
                t =0;
       
            case '2'
                disp ('Right click to reset the ONSET of the Warning Signal \n');
                button = 1;
                while sum(button) <=1   % read ginputs until a mouse right-button occurs
                   [x,y,button] = ginput(1);
                end

                plot(x,y, 'm*');
                new_onset  = round(x);
                fprintf('old onset: %d, new onset: %d \n', firstSig_indx, new_onset);
                warning_onset = new_onset;
                t =0;
                if ~isempty(warning_marker_index)
                    warning_marker_index(n) = buffered_index+warning_onset;
                else
                    warning_marker_index = buffered_index+warning_onset;
                end
                
            otherwise 
                prompt = 'Wrong input. 2:yes 4:no \n';
                answer = input(prompt,'s');

        end
    end
    n=n+1;
    close all
end
    close all
    
    
    
%     %for debugging purposes
%     %this plot shows the markers together with the peaks detected for each
%     %signal. Written for just one index marker
%     figure (2)
%     plot(EMG_data_marker);
%     hold on;
%     %xline(index_markers, 'r')
%     for i=1:length(index_markers)
%     plot([index_markers(i),index_markers(i)],[-max(EMG_data_marker),max(EMG_data_marker)],'r');
%     hold on;
%     end
%     hold on
%     plot(index+loc, EMG_data_marker(index+loc),'*g');



%close all
fprintf('Number of Trials: %d \n',length(warning_marker_index));

end