function [Epochs_pre, params, trigger_order, run_loc, task_pre] = computeModel(subject)

%%%%%%%%%%%%%%%%%%%%
%% Initialization %%
%%%%%%%%%%%%%%%%%%%%
clearvars -except subject;
 clc; rng('default');


%%%%%%%%%%%%%%%%%%
%% Load Dataset %%
%%%%%%%%%%%%%%%%%%
currentPath = pwd;
dataPath = [currentPath '/data'];
[task_pre, run_loc] = loadData_v2(subject, ['data/subject' num2str(subject) '/*.gdf']);
params = setParameters(task_pre(1).header);
params.triggerChannel = size(task_pre(1).data, 2);
params.eegChannels_noFCz = [1:64];
channels = length(params.eegChannels_noFCz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocess EEG signals %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%temporal filter
[b,a] = butter(params.spectralFilter_theta.order, params.spectralFilter_theta.freqs./(params.fsamp/2), 'bandpass');
params.spectralFilter.a = a;
params.spectralFilter.b = b;
task_pre.data(:, [params.eegChannels_noFCz]) = filtfilt(b, a, task_pre.data(:, [params.eegChannels_noFCz]));

%%%%%%%%%%%%%%
%% Epoching (pre stim) %%
%%%%%%%%%%%%%%  
trigger = task_pre.data(:,end);
if (isempty(task_pre.data) == 0)
    %Number of car movements
    [Positive_triggers, Negative_triggers, Neutral_triggers] = computeIndex(trigger);
    trigger_order = [Positive_triggers, Negative_triggers, Neutral_triggers];
    n_training_trials = length(Positive_triggers) + length(Negative_triggers) + length(Neutral_triggers);

    Epochs_pre = struct('data', nan(length(params.epochSample), length(params.eegChannels_noFCz), n_training_trials),...
        'label', nan(n_training_trials, 1),...
        'fileID', nan(n_training_trials, 1),...
        'sessionID', nan(n_training_trials, 1));

        for i_trial = 1:length(Positive_triggers)
            Epochs_pre.data(:, :, i_trial) = task_pre.data(Positive_triggers(i_trial) + params.epochSample, params.eegChannels_noFCz);
            temp = find(Positive_triggers(i_trial) <= task_pre.eof_eeg, 1, 'first');
            Epochs_pre.fileID(i_trial) = temp;
            Epochs_pre.label(i_trial) = 1;
        end
        for i_trial = 1:length(Negative_triggers)
            Epochs_pre.data(:, :, length(Positive_triggers)+i_trial) = task_pre.data(Negative_triggers(i_trial) + params.epochSample, params.eegChannels_noFCz);
            temp = find(Negative_triggers(i_trial) <= task_pre.eof_eeg, 1, 'first');
            Epochs_pre.fileID(length(Positive_triggers)+i_trial) = temp;
            Epochs_pre.label(length(Positive_triggers)+i_trial) = 2;
        end
        for i_trial = 1:length(Neutral_triggers)
            Epochs_pre.data(:, :, length(Positive_triggers)+length(Negative_triggers)+i_trial) = task_pre.data(Neutral_triggers(i_trial) + params.epochSample, params.eegChannels_noFCz);
            temp = find(Neutral_triggers(i_trial) <= task_pre.eof_eeg, 1, 'first');
            Epochs_pre.fileID(length(Positive_triggers)+length(Negative_triggers)+i_trial) = temp;
            Epochs_pre.label(length(Positive_triggers)+length(Negative_triggers)+i_trial) = 3;
        end
    end
end

function [Positive_triggers, Negative_triggers, Neutral_triggers] = computeIndex(trigger)
    %loop through car movement triggers, if followed by 102 or 103, then remove
    %it
    Positive_triggers = [];
    Negative_triggers = [];
    Neutral_triggers = [];
    trigger_positions = 1:length(trigger);
    delete_index = [find(trigger == 0); find(trigger == 1)];
    trigger_positions(delete_index) = [];
    trigger(delete_index) = [];
    for i = 1:length(trigger)-1
        if trigger(i) == 100 && trigger(i+1) == 102
            Positive_triggers = cat(2,Positive_triggers, trigger_positions(i));
        elseif trigger(i) == 100 && trigger(i+1) == 103
            Negative_triggers = cat(2,Negative_triggers, trigger_positions(i));
        elseif trigger(i) == 100 && trigger(i+1) == 100
            Neutral_triggers = cat(2,Neutral_triggers, trigger_positions(i));
        end
    
    end
end
