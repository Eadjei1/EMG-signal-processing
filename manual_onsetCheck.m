function [scm_EMGs,scm_emg_threshold,SCM_Time_output,Time_output]=manual_onsetCheck(headerFile,filename,outputData,outputData_rect,marker_index, ...
    onset_duration,maxDuration,order,sound,emg_prelength, scm_preLength,fs)
%%
%initializing and setting up the stage
ind = find (headerFile == '.');
outputName = strcat(headerFile(1:ind-1),'_',sound, '_res.mat');
% soundWave = outputData(1,:);
soundWave_rect = outputData_rect(1,:);
EMGs = outputData(2:end,:);
EMGs_rect = outputData_rect(2:end,:);

%finding the primary and secondary muscles for Open and Lift
pat1 = "lift";
pat2 = "open";
prim_Open = order == "EDC";
prim_Lift = order == "IDL";
sec1_mus = order == "BIC";
sec2_mus = order == "TRI";
wristexten = order =="ECR";
wristflex = order =="FCR";
fingerflex = order =="FDS";
otherBIC = order =="OtherBIC";

if any(contains(order, 'SCM')) 
    scmr_mus = order == "SCMr";
    scml_mus = order == "SCMl";
    SCMr = EMGs(scmr_mus,:);
    SCMl = EMGs(scml_mus,:);
    SCMr_rect = EMGs_rect(scmr_mus,:);
    SCMl_rect = EMGs_rect(scml_mus,:);
end

BIC = EMGs(sec1_mus,:);
TRI = EMGs(sec2_mus,:);
IDL = EMGs(prim_Lift,:);
EDC = EMGs(prim_Open,:);
ECR = EMGs(wristexten,:);
FCR = EMGs(wristflex,:);
FDS = EMGs(fingerflex,:);
OtherBIC = EMGs(otherBIC,:);

BIC_rect = EMGs_rect(sec1_mus,:);
TRI_rect = EMGs_rect(sec2_mus,:);
IDL_rect = EMGs_rect(prim_Lift,:);
EDC_rect = EMGs_rect(prim_Open,:);
ECR_rect = EMGs_rect(wristexten,:);
FCR_rect = EMGs_rect(wristflex,:);
FDS_rect = EMGs_rect(fingerflex,:);
OtherBIC_rect = EMGs_rect(otherBIC,:);

if contains(filename, pat2,"IgnoreCase",true)
    emg_muscle = find(prim_Open == 1);
    %Lift+Open
    if contains(filename,pat1, "IgnoreCase",true)
        emg_muscle = [find(prim_Lift == 1),find(prim_Open == 1)];

    end
%for lift only
elseif contains(filename,pat1, "IgnoreCase",true)
    emg_muscle = find(prim_Lift == 1);
end

Time_output = zeros(length(marker_index),length(emg_muscle)*2);
Time_output_ms = zeros(length(marker_index),length(emg_muscle)*2);
Actual_onoff = zeros(length(marker_index),length(emg_muscle)*2);
SI_on_off = zeros(length(marker_index),length(emg_muscle)*2);
SCM_Time_output = zeros(length(marker_index),4); %changed column to 4 instead of 2
ActualSCM_onoff = zeros(length(marker_index),4);
scm_EMGs = zeros(length(marker_index),2);
scm_emg_threshold = zeros(length(marker_index),1);
ECR_output = zeros(length(marker_index),2);
% ECR_RMS = zeros(length(marker_index),size(EMGs,1));
% ECR_RMS1 =  zeros(length(marker_index),size(EMGs,1));

% ECR_onset = 0;
% ECR_onset_sec = 0;
%%
channel = 1;
coef = 5;
SCM_coef = 3;
chan = 1;
prelength = 0.04*fs; %20ms after the sound cue
if emg_prelength == 0
    disp("this is CCC")
    prelength = 2*fs;
end

while channel <= length(emg_muscle)+1
    prim_muscle = emg_muscle(chan);
    %plotting either IDL or EDC as second row
    if order(prim_muscle) == "IDL"
        muscle = EDC; 
        muscle_rect = EDC_rect;
        Label = order(prim_Open);
        disp('primary muscle is the IDL')
    else
        muscle = IDL;
        muscle_rect = IDL_rect;
        Label = order(prim_Lift);
        disp('primary muscle is the EDC')
    end

    trials = 1;
    num = emg_muscle(chan);
    signal = EMGs(num,:);
    signal_rect = EMGs_rect(num,:);
    for marker = marker_index
                %finding threshold of onset detection using the envelope
                mean_emg = mean(signal_rect(marker-emg_prelength:marker+prelength));
                std_emg = std(signal_rect(marker-emg_prelength:marker+prelength));
                threshold = mean_emg+(coef*std_emg);
                
                if (marker+onset_duration+maxDuration) <= length(signal_rect)
                    emg_interval = marker-onset_duration:marker+onset_duration+maxDuration;
                    emg_interval_sec = linspace(0-(onset_duration/fs), (onset_duration+maxDuration)/fs, length(emg_interval));
                else
                    emg_interval = marker:length(signal);
                    emg_interval_sec = linspace(0-(onset_duration/fs), (onset_duration+maxDuration)/fs, length(emg_interval));
                    disp('the end has come')
                end
                
%                 %marker2onset_interval = marker:marker+onset_duration+maxDuration;
%                 above_threshold = abs(signal_rect(emg_interval)) > threshold;
%                 %emg_loc_4rm_marker = find(above_threshold);
%                 emg_loc = find(above_threshold);
%                 %emg_loc = emg_loc_4rm_marker+onset_duration;
% 
%                 % Find consecutive regions where the logical array is true
% %                 regions = regionprops(above_threshold, 'Area', 'PixelIdxList');
% %                 above_threshold_regions = [regions.Area];
% %                 above_threshold_indices = {regions.PixelIdxList};
%                 
%                 duration_threshold = 0.02 * fs; % 10ms in terms of samples (assuming sampling frequency is 'fs')
%                 
%                 onset_index = 0; % Initialize the onset index to zero
% %                 offset_index = 0; % Initialize the onset index to zero
%                 % Iterate through the regions and find the first region that exceeds the duration threshold
%                 for i = 1:numel(above_threshold_regions)
%                     if above_threshold_regions(i) >= duration_threshold
%                         onset_index = above_threshold_indices{i}(1)+onset_duration; % Set the onset index to the first index of the region
%                         %offset_index = above_threshold_indices{i}(end); % Set the onset index to the first index of the region
%                        
%                         break; % Exit the loop as we have found the onset
%                     end
%                 end
%                  for i = 1:length(emg_loc)
%                     start_index = emg_loc(i);
%                     end_index = start_index + duration_threshold - 1;
%                     
%                     if end_index <= length(signal_rect(emg_interval)) && all(abs(signal_rect(start_index:end_index)) > threshold)
%                         onset_index = start_index;
%                         break; % Exit the loop as we have found the onset
%                         disp("should not")
%                     end
%                 end
                
                emg_loc = find(abs(signal_rect(emg_interval)) > threshold);
    
                onset = marker-onset_duration+emg_loc(1);
                onset_sec = (emg_loc(1)-onset_duration)/fs;
                offset = marker-onset_duration+emg_loc(end);
                offset_sec = (emg_loc(end)-onset_duration)/fs;
%                 Sound_cur = soundWave(emg_interval);
                Sound_cur_rect = soundWave_rect(emg_interval);
                
                figure
                subplot(5, 1, 1)
%                 plot(emg_interval_sec,signal(emg_interval))
%                 hold on
                plot(emg_interval_sec,signal_rect(emg_interval))
                title(strcat(filename,' ', sound, ': Trial Number: ',num2str(trials)));
    
                hold on
                plot(onset_sec, signal(onset),'*g');
                plot(offset_sec, signal(offset),'*k');
                plot([0,0],[-max(signal(emg_interval)),max(signal(emg_interval))],'r');
                legend(order(prim_muscle),'onset', 'offset', sound)
  
    
                subplot(5, 1, 2)
%                 plot(emg_interval_sec,muscle(emg_interval))
%                 hold on
                plot(emg_interval_sec,muscle_rect(emg_interval))
                hold on
                plot([0,0],[-max(muscle(emg_interval)),max(muscle(emg_interval))],'r');
                legend(Label,sound)
                
                subplot(5,1,3)
                plot(emg_interval_sec,BIC_rect(emg_interval))
                hold on
%                 plot(emg_interval_sec,BIC_rect(emg_interval),'c')
                plot([0,0],[-max(BIC(emg_interval)),max(BIC(emg_interval))],'r');
                legend('BIC','rectified',sound)
                
                subplot(5, 1, 4)
                plot(emg_interval_sec,TRI_rect(emg_interval))
                hold on
%                 plot(emg_interval_sec,TRI_rect(emg_interval),'c')
                plot([0,0],[-max(TRI(emg_interval)),max(TRI(emg_interval))],'r');
                legend('TRI','rectified',sound)
    
                subplot(5, 1, 5)
                plot(emg_interval_sec,Sound_cur_rect)
%                 hold on
%                 plot(emg_interval_sec,Sound_cur_rect,'c')
                legend('Sound Signal', 'rectified')
                
                prompt = 'Do you want reset the onset or offset or both? 1:both, 2:onset, 3:offset, 4: none \n';
                answer = input(prompt,'s');
                t = 1;
                while t ==1
                    switch answer 
                        case '4'
                            disp ('You have decided not to reset onset and offset, proceed to the next');
                            t =0;
                        case '1'
                            disp ('Right click to reset the ONSET of EMG \n');
                            button = 1;
                            while sum(button) <=1   % read ginputs until a mouse right-button occurs
                               [x,y,button] = ginput(1);
                            end
    
                            plot(x,y, 'm*');
                            %new_onset  = (x*fs)+marker;
                            new_onset  = (x*fs)+marker;
                            new_onset_sec = x;
                            fprintf('old onset: %d, new onset: %d \n', onset, new_onset);
                            fprintf('old onset_sec: %d, new onset_sec: %d', onset_sec, new_onset_sec);
                            onset = new_onset;
                            onset_sec = new_onset_sec;
    
                            disp ('Right click to reset the OFFSET of EMG');
                            button = 1;
                            while sum(button) <=1   % read ginputs until a mouse right-button occurs
                               [x,y,button] = ginput(1);
                            end
                            plot(x,y, 'm*');
                            offset  = (x*fs)+marker;
                            offset_sec = x;
                            %offset  = marker-(1*fs)+x;
                            t =0;
                        case '2'
                            disp ('Right click to reset the ONSET of EMG \n');
                            button = 1;
                            while sum(button) <=1   % read ginputs until a mouse right-button occurs
                               [x,y,button] = ginput(1);
                            end
    
                            plot(x,y, 'm*');
                            %new_onset  = marker-(1*fs)+x;
%                             new_onset  = (x*fs)+marker;
%                             fprintf('old onset: %d, new onset: %d', onset, new_onset);
%                             onset = new_onset;
                            new_onset  = (x*fs)+marker;
                            new_onset_sec = x;
                            fprintf('old onset: %d, new onset: %d/n', onset, new_onset);
                            fprintf('old onset_sec: %d, new onset_sec: %d', onset_sec, new_onset_sec);
                            onset = new_onset;
                            onset_sec = new_onset_sec;
                            t =0;
                    
                        case '3'
                            disp ('Right click to reset the OFFSET of EMG \n');
                            button = 1;
                            while sum(button) <=1   % read ginputs until a mouse right-button occurs
                               [x,y,button] = ginput(1);
                            end
                            plot(x,y, 'm*');
                            offset  = (x*fs)+marker;
                            offset_sec = x;
                            %offset  = marker-(1*fs)+x;
                            t =0;
                        otherwise 
                            prompt = 'Wrong input. 1:both, 2:onset, 3:offset, 4: none \n';
                            answer = input(prompt,'s');
                            
                    end
                end
    
                close all
                disp ('Yaay, you are done selecting the IDL(EDC) onset and offset for this marker!');

                %% adding ecr
                %checking for ECR onsets as well since hand opening
                %encompasses both
                if channel == 1 % just so you check the ECR only once
%                     mean_ecr_emg = mean(ECR(marker:marker+scm_preLength));
%                     std_ecr_emg = std(ECR(marker:marker+scm_preLength));
%                     ecr_emg_threshold(trials,channel) = mean_ecr_emg+(coef*std_ecr_emg);
%                   
                    subplot(5, 1, 1)
                    plot(emg_interval_sec,ECR_rect(emg_interval))
                    title (strcat('ECR onset '' ', filename,' ',sound,': Trial Number: ',num2str(trials)));
            
                    hold on
%                     plot(emg_interval_sec,ECR_rect(emg_interval),'c')
                    plot([0,0],[-max(ECR(emg_interval)),max(ECR(emg_interval))],'r');
                    legend('ECR',sound)
        
                    subplot(5, 1, 2)
                    plot(emg_interval_sec,EDC_rect(emg_interval))
                    hold on
%                     plot(emg_interval_sec,EDC_rect(emg_interval),'c')
                    plot([0,0],[-max(EDC(emg_interval)),max(EDC(emg_interval))],'r');
                    legend('EDC',sound)

                    subplot(5, 1, 3)
                    plot(emg_interval_sec,FDS_rect(emg_interval))
                    hold on
%                     plot(emg_interval_sec,FDS_rect(emg_interval),'c')
                    plot([0,0],[-max(FDS(emg_interval)),max(FDS(emg_interval))],'r');
                    legend('FDS',sound)

                    subplot(5, 1, 4)
                    plot(emg_interval_sec,FCR_rect(emg_interval))
                    hold on
%                     plot(emg_interval_sec,FCR_rect(emg_interval),'c')
                    plot([0,0],[-max(FCR(emg_interval)),max(FCR(emg_interval))],'r');
                    legend('FCR',sound)
        
%                     subplot(5, 1, 5)
%                     plot(emg_interval_sec,Sound_cur_rect)
% %                     hold on
% %                     plot(emg_interval_sec,Sound_cur_rect,'c')
%                     legend('Sound Signal')
                    subplot(5, 1, 5)
                    plot(emg_interval_sec,OtherBIC_rect(emg_interval))
                    hold on
%                     plot(emg_interval_sec,FCR_rect(emg_interval),'c')
                    plot([0,0],[-max(OtherBIC_rect(emg_interval)),max(OtherBIC_rect(emg_interval))],'r');
                    legend('OtherBIC',sound)
                    
                    prompt = 'Do you want reset the ECR onset or offset or both? 2:onset, 4: none \n';
                    answer = input(prompt,'s');
                    t = 1;
                    while t ==1
                        switch answer 
                            case '4'
                                disp ('You have decided not to reset onset and offset, proceed to the next');
                                ECR_onset = 0;
                                ECR_onset_sec = 0;
                                t =0;
                            case '2'
                                disp ('Right click to reset the ONSET of EMG \n');
                                button = 1;
                                while sum(button) <=1   % read ginputs until a mouse right-button occurs
                                   [x,y,button] = ginput(1);
                                end
        
                                plot(x,y, 'm*');
                                new_ECRonset  = (x*fs)+marker;
                                fprintf('new onset: %d', new_ECRonset);
                                ECR_onset = new_ECRonset;
                                ECR_onset_sec = x;
                                t =0;
                            otherwise 
                                prompt = 'Wrong input. 1:both, 2:onset, 3:offset, 4: none \n';
                                answer = input(prompt,'s');
                                
                        end
                    end
                    close all
                else
                    disp('You have already selected ECR on and off for this trial')
                end
                
                %end of addition
                %% 
                time_onset = ((onset-marker)/fs)*1000; %converting to ms
                time_offset = ((offset-marker)/fs)*1000;
                time_onset_ms = onset_sec*1000;
                time_offset_ms = offset_sec*1000;

                time_ECRonset = ((ECR_onset-marker)/fs)*1000; %converting to ms
                time_ECRonset_ms = ECR_onset_sec*1000;

                actual_time_onset = (onset/fs); %converting to seconds
                actual_time_offset = (offset/fs);
                
                Time_output(trials,channel:channel+1) = [time_onset, time_offset];
                Time_output_ms(trials,channel:channel+1) = [time_onset_ms, time_offset_ms];
                ECR_output(trials,channel:channel+1) = time_ECRonset;
                ECR_output_ms(trials,channel:channel+1) = time_ECRonset_ms;
                SI_on_off(trials,channel:channel+1) = [actual_time_onset, actual_time_offset];
                Actual_onoff(trials,channel:channel+1) = [onset, offset];
    
                close all

                %checking for SCM activation 
                if any(contains(order, 'SCM'))
                    if channel == 1 % just so you check the SCM only once
                        mean_scm_emg = mean(SCMl(marker-scm_preLength:marker));
                        std_scm_emg = std(SCMl(marker-scm_preLength:marker));
                        scm_emg_threshold(trials,channel) = mean_scm_emg+(SCM_coef*std_scm_emg);
        
                        subplot(4, 1, 1)
                        plot(emg_interval_sec,SCMr_rect(emg_interval))
                        title (strcat('SCM+ Check for '' ', filename,' ',sound,': Trial Number: ',num2str(trials)));
                
                        hold on
%                         plot(emg_interval_sec,SCMr_rect(emg_interval),'c')
                        plot([0,0],[-max(SCMr(emg_interval)),max(SCMr(emg_interval))],'r');
                        legend('SCMr',sound)
            
                        subplot(4, 1, 2)
                        plot(emg_interval_sec,SCMl_rect(emg_interval))
                        hold on
%                         plot(emg_interval_sec,SCMl_rect(emg_interval),'c')
                        plot([0,0],[-max(SCMl(emg_interval)),max(SCMl(emg_interval))],'r');
                        legend('SCMl',sound)
            
                        subplot(4, 1, 3)
%                         plot(emg_interval_sec,Sound_cur_rect)
% %                         hold on
% %                         plot(emg_interval_sec,Sound_cur_rect,'c')
%                         legend('Sound Signal')
                        
                        plot(emg_interval_sec,OtherBIC_rect(emg_interval))
                        hold on
    %                     plot(emg_interval_sec,FCR_rect(emg_interval),'c')
                        plot([0,0],[-max(OtherBIC_rect(emg_interval)),max(OtherBIC_rect(emg_interval))],'r');
                        legend('OtherBIC',sound)
                        
                        subplot(4, 1, 4)
                        plot(emg_interval_sec,Sound_cur_rect)
        %                 hold on
        %                 plot(emg_interval_sec,Sound_cur_rect,'c')
                        legend('Sound Signal', 'rectified')

                        prompt = 'Select SCM on/off 1:Both 2:Right 3: Left  4: No SCM+ \n';
                        answer = input(prompt,'s');
                        t = 1;
                        while t ==1
                            switch answer 
                                case '4'
                                    t =0;
                                case '1'
                                    disp ('Right click to reset the ONSET of right SCM \n');
                                    button = 1;
                                    while sum(button) <=1   % read ginputs until a mouse right-button occurs
                                       [x,y,button] = ginput(1);
                                    end
            
                                    plot(x,y, 'm*');
                                    SCMr_onset  = (x*fs)+marker;
                                    
                                    disp ('Right click to reset the OFFSET of right SCM');
                                    button = 1;
                                    while sum(button) <=1   % read ginputs until a mouse right-button occurs
                                       [x,y,button] = ginput(1);
                                    end
                                    plot(x,y, 'm*');
                                    SCMr_offset  = (x*fs)+marker;
                                    t =0;

                                    disp ('Right click to reset the ONSET of left SCM \n');
                                    button = 1;
                                    while sum(button) <=1   % read ginputs until a mouse right-button occurs
                                       [x,y,button] = ginput(1);
                                    end
            
                                    plot(x,y, 'm*');
                                    SCMl_onset  = (x*fs)+marker;
                                    
                                    disp ('Right click to reset the OFFSET of left SCM');
                                    button = 1;
                                    while sum(button) <=1   % read ginputs until a mouse right-button occurs
                                       [x,y,button] = ginput(1);
                                    end
                                    plot(x,y, 'm*');
                                    SCMl_offset  = (x*fs)+marker;
                                    t =0;

                                case '2'
                                    disp ('Right click to reset the ONSET of right SCM \n');
                                    button = 1;
                                    while sum(button) <=1   % read ginputs until a mouse right-button occurs
                                       [x,y,button] = ginput(1);
                                    end
            
                                    plot(x,y, 'm*');
                                    SCMr_onset  = (x*fs)+marker;
                                    
                                    disp ('Right click to reset the OFFSET of right SCM');
                                    button = 1;
                                    while sum(button) <=1   % read ginputs until a mouse right-button occurs
                                       [x,y,button] = ginput(1);
                                    end
                                    plot(x,y, 'm*');
                                    SCMr_offset  = (x*fs)+marker;
                                    t =0;

                                    SCMl_offset = 0; SCMl_onset = SCMl_offset;

                                case '3'

                                    disp ('Right click to reset the ONSET of left SCM \n');
                                    button = 1;
                                    while sum(button) <=1   % read ginputs until a mouse right-button occurs
                                       [x,y,button] = ginput(1);
                                    end
            
                                    plot(x,y, 'm*');
                                    SCMl_onset  = (x*fs)+marker;
                                    
                                    disp ('Right click to reset the OFFSET of left SCM');
                                    button = 1;
                                    while sum(button) <=1   % read ginputs until a mouse right-button occurs
                                       [x,y,button] = ginput(1);
                                    end
                                    plot(x,y, 'm*');
                                    SCMl_offset  = (x*fs)+marker;
                                    t =0;

                                    SCMr_offset = 0; SCMr_onset = SCMr_offset;
                                   
                                
                                otherwise 
                                    prompt = 'Wrong input. Is there SCM+? Select on/offset? 1:Yes, both, 2:Yes,onset 4: No SCM+ \n';
                                    answer = input(prompt,'s');
                                    
                            end
                        end
                        close all
                    else
                        disp('You have already selected SCM on and off for this trial')
                    end
                end

                if exist('SCMr_onset','var') || exist('SCMl_onset','var')
                    SCMr_time_onset = ((SCMr_onset-marker)/fs)*1000; %converting to ms
                    SCMr_time_offset = ((SCMr_offset-marker)/fs)*1000;

                    SCMl_time_onset = ((SCMl_onset-marker)/fs)*1000; %converting to ms
                    SCMl_time_offset = ((SCMl_offset-marker)/fs)*1000;

                    SCM_Time_output(trials,:) = [SCMr_time_onset, SCMr_time_offset, SCMl_time_onset, SCMl_time_offset];
                    ActualSCM_onoff(trials,:) = [SCMr_onset,SCMr_offset,SCMl_onset, SCMl_offset];
                    scm_EMGs(trials,channel:channel+1) = max(abs(EMGs(end-1:end,SCMl_onset:SCMl_offset)), [],2)';
                else
                    disp('No SCM activity present. Use onset of primary muscle')
                end

                if ~any(contains(order, 'SCM'))
                     save(outputName, 'Time_output', 'Time_output_ms','Actual_onoff', 'SI_on_off')
                else
                    %save(outputName, 'Time_output','Time_output_ms','Actual_onoff' ,'SI_on_off', ...
                    %    'SCM_Time_output','ActualSCM_onoff','scm_emg_threshold', 'scm_EMGs')

                    save(outputName, 'Time_output','Time_output_ms','Actual_onoff' ,'SI_on_off', ...
                        'SCM_Time_output','ActualSCM_onoff','scm_emg_threshold','scm_EMGs')
                end
                clear SCMl_onset SCMr_onset
                trials=trials+1;
    end
    channel = channel +2;
    chan = chan+1;

end
  
end