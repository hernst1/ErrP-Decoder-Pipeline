
%% (1) Trials extraction and subject wise plotting
close all
clear
addpath(genpath('../PL_enhance_BCI/code/functions'))
subjectID = ['4' '6']; %['1' '2' '3' '5'];% 
amplitude_threshold = 10000;
load('chanlocs64.mat')
channellocations = chanlocs;
channellocations(32) = [];
chan = [37, 46, 47]; %[5, 39];

deep_colors = {[31,120,180], [227,26,28], [51,160,44]};
deep_colors = cellfun(@(x) x./255, deep_colors, 'UniformOutput', false);

currentPath = pwd;
dataPath = [currentPath '/data'];

eegEpochs = struct;
triggers = struct;
run_locs = cell(length(subjectID), 1);

for subID = 1:length(subjectID)
        [eegEpochs_pre, params, trigger_order, run_loc, task_pre] = computeModel(subjectID(subID));

        trigger_unorder = [trigger_order' eegEpochs_pre.label];
        trigger_order = sortrows([trigger_order' eegEpochs_pre.label],1);
        propagate_run = zeros(size(trigger_order,1),1);
        
        for run = length(run_loc):-1:1
            indices = find(trigger_order(:,1) <= run_loc(run));
            propagate_run(indices) = run;
        end
        
        trigger_order = [propagate_run trigger_order];
        triggers.subjectID(subID).ordered = trigger_order;
        triggers.subjectID(subID).unordered = trigger_unorder;

        % filename = fullfile(dataPath, ['trigger_order_sub' num2str(subjectID(subID)) '.mat']);
        % save(filename, 'trigger_order');
        % disp(['Saved: ', filename]);

        run_locs{subID} = run_loc;
        eegEpochs_pre.data(:,32,:) = [];
        %channellocations(32) = [];
        signal_range = dsearchn(params.epochTime', 0.0):dsearchn(params.epochTime', 1.0);
        
        for chan_id = 1:length(chan)
            %figure;
            lineLegend = [];

            % UNCOMMENT for mu power analysis and baseline normalization
            % baseline_range = dsearchn(params.epochTime', -0.5):dsearchn(params.epochTime', 0.0);
            % mu_power = (eegEpochs_pre.data).^2;
            % baseline_mu = mean(mu_power(baseline_range, :, :),1);
            %

            signal_range = dsearchn(params.epochTime', 0.0):dsearchn(params.epochTime', 1.0);
            keep_index = not(squeeze(any(any(abs(eegEpochs_pre.data(signal_range,:,:)) > amplitude_threshold))));
            
            logical_index = eegEpochs_pre.label == 1;
            logical_index_positive = logical_index & keep_index; disp(['total positive valence trials:' string(sum(logical_index_positive))])
            logical_index = eegEpochs_pre.label == 2;
            logical_index_negative = logical_index & keep_index; disp(['total negative valence trials:' string(sum(logical_index_negative))])
            logical_index = eegEpochs_pre.label == 3;
            logical_index_neutral = logical_index & keep_index; disp(['total neutral trials:' string(sum(logical_index_neutral))])
            
            %temp = eegEpochs_pre.data;
            %temp = (eegEpochs_pre.data).^2; % mu power
            
            % Plotting negative class %
            % Theta %
            %logical_index_negative = eegEpochs_pre.label == 2 & keep_index;
            %shadedLine2 = shadedErrorBar(params.epochTime, nanmean(temp(:, chan(chan_id), logical_index_negative), 3), nanstd(temp(:, chan(chan_id), logical_index_negative), [], 3)/sqrt(length(logical_index_negative)), {'Color', deep_colors{1}, 'Linewidth', 3, 'DisplayName', 'Negative valence', 'lineWidth', 3}, 0.1);hold on
            %plot(params.epochTime, nanmean(temp(:,chan(chan_id), logical_index_negative), 3), 'DisplayName', 'Negative valence trials', 'lineWidth', 3); hold on
            
            % Mu baseline normalization %
            % baseline_mu_negative = mean(baseline_mu(:, chan(chan_id), logical_index_negative),3);
            % norm_temp = (temp(:, chan(chan_id), logical_index_negative) - baseline_mu_negative)/baseline_mu_negative;
            % shadedLine2 = shadedErrorBar(params.epochTime, nanmean(norm_temp, 3), nanstd(norm_temp, [], 3)/sqrt(length(logical_index_negative)), {'Color', deep_colors{1}, 'Linewidth', 3, 'DisplayName', 'Negative valence', 'lineWidth', 3}, 0.1);hold on
            
            %lineLegend = cat(1, lineLegend, shadedLine2.mainLine);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Plotting positive class %
            % Theta %
            %logical_index_positive = eegEpochs_pre.label == 1 & keep_index;
            %shadedLine1 = shadedErrorBar(params.epochTime, nanmean(temp(:, chan(chan_id), logical_index_positive), 3), nanstd(temp(:, chan(chan_id), logical_index_positive), [], 3)/sqrt(length(logical_index_positive)), {'Color', deep_colors{2}, 'Linewidth', 3, 'DisplayName', 'Positive valence', 'lineWidth', 3}, 0.1);hold on
            %plot(params.epochTime, nanmean(temp(:,chan(chan_id), logical_index_positive), 3), 'DisplayName', 'Positive valence trials', 'lineWidth', 3); hold on
            
            % Mu baseline normalization %
            % baseline_mu_positive = mean(baseline_mu(:, chan(chan_id), logical_index_positive),3);
            % norm_temp = (temp(:, chan(chan_id), logical_index_positive) - baseline_mu_positive)/baseline_mu_positive;
            % shadedLine1 = shadedErrorBar(params.epochTime, nanmean(norm_temp, 3), nanstd(norm_temp, [], 3)/sqrt(length(logical_index_positive)), {'Color', deep_colors{2}, 'Linewidth', 3, 'DisplayName', 'Positive valence', 'lineWidth', 3}, 0.1);hold on
            
            %lineLegend = cat(1, lineLegend, shadedLine1.mainLine);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Plotting neutral class %
            % Theta %
            %logical_index_neutral = eegEpochs_pre.label == 3 & keep_index;
            %shadedLine3 = shadedErrorBar(params.epochTime, nanmean(temp(:, chan(chan_id), logical_index_neutral), 3), nanstd(temp(:, chan(chan_id), logical_index_neutral), [], 3)/sqrt(length(logical_index_neutral)), {'Color', deep_colors{3}, 'Linewidth', 3, 'DisplayName', 'Neutral', 'lineWidth', 3}, 0.1);hold on
            %plot(params.epochTime, nanmean(temp(:,chan(chan_id), logical_index_neutral), 3), 'DisplayName', 'Neutral trials', 'lineWidth', 3); hold on
            
            % Mu baseline normalization %
            % baseline_mu_neutral = mean(baseline_mu(:, chan(chan_id), logical_index_neutral),3);
            % norm_temp = (temp(:, chan(chan_id), logical_index_neutral) - baseline_mu_neutral)/baseline_mu_neutral;
            % shadedLine3 = shadedErrorBar(params.epochTime, nanmean(norm_temp, 3), nanstd(norm_temp, [], 3)/sqrt(length(logical_index_neutral)), {'Color', deep_colors{3}, 'Linewidth', 3, 'DisplayName', 'Neutral', 'lineWidth', 3}, 0.1);hold on

            %lineLegend = cat(1, lineLegend, shadedLine3.mainLine);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % title(['ErrPs at ' channellocations(chan(chan_id)).labels ' average']); grid on; grid minor;
            % %title(['Mu Power at ' channellocations(chan(chan_id)).labels ' average']); grid on; grid minor;
            % legend(lineLegend);
            % 
            % if chan_id == 2
            %     axis([-0.5 1.5 -4 4])%max(max(nanmean(temp(:, chan(2), logical_index_positive), 3), nanmean(temp(:, chan(2), logical_index_negative), 3)))]);
            % else
            %     axis([-0.5 1.5 -6 7])%max(max(nanmean(temp(:, chan(1), logical_index_positive), 3), nanmean(temp(:, chan(1), logical_index_negative), 3)))]);
            % end
            % xlabel('Time [s]', 'FontSize',10);
            % ylabel('Amplitude [\muV]', 'FontSize',10);
            % %ylabel('Power [\muV]^2', 'FontSize',10);
            % h = plot([0 0], get(gca, 'YLim'), '--k', 'lineWidth', 3.0);
            % h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % set(gca, 'FontSize', 12);
        end
        
        % Select time-points and 'removed' ch for topoplots
        % figure;
        % time_index = [0.35, 0.5, 0.70];
        % color_ticks = [-3 0 3];
        % r = 1;
        % eegEpochs_temp = eegEpochs_pre;
        % if strcmp(subjectID, '6') == 1
        %     eegEpochs_temp.data(:,[13,19],:) = [];
        %     channellocations([13,19]) = [];
        % end

        % Mu Band Topoplots
        % exEEG = mean(eegEpochs_temp.data(:,:,logical_index_negative).^2, 3);
        % subplot(1,3,r)
        % time_samples = dsearchn(params.epochTime', params.epochTime(end));
        % topoplot(exEEG(time_samples,:), channellocations, 'MapLimits', [min(color_ticks) max(color_ticks)]);
        % title(['Negative: 0 - ' num2str(params.epochTime(end), '%.2f') ' s']);
        % colorbar('Ticks',min(color_ticks):5:max(color_ticks));
        % r = r + 1;
        % 
        % exEEG = mean(eegEpochs_temp.data(:,:,logical_index_positive).^2, 3);
        % subplot(1,3,r)
        % time_samples = dsearchn(params.epochTime', params.epochTime(end));
        % topoplot(exEEG(time_samples,:), channellocations, 'MapLimits', [min(color_ticks) max(color_ticks)]);
        % title(['Positive: 0 - ' num2str(params.epochTime(end), '%.2f') ' s']);
        % r = r + 1;
        % colorbar('Ticks',min(color_ticks):5:max(color_ticks));
        % 
        % exEEG = mean(eegEpochs_temp.data(:,:,logical_index_neutral).^2, 3);
        % subplot(1,3,r)
        % time_samples = dsearchn(params.epochTime', params.epochTime(end));
        % topoplot(exEEG(time_samples,:), channellocations, 'MapLimits', [min(color_ticks) max(color_ticks)]);
        % title(['Neutral: 0 - ' num2str(params.epochTime(end), '%.2f') ' s']);
        % r = r + 1;
        % colorbar('Ticks',min(color_ticks):5:max(color_ticks));
        % 
        % sgtitle("Inspecting Frontal Asymmetry in Mu")
        
        % Theta Band Topoplots
        % exEEG = mean(eegEpochs_temp.data(:,:,logical_index_negative), 3);
        % for i_time = 1:length(time_index)
        %     subplot(3,length(time_index),r)
        %     time_samples = dsearchn(params.epochTime', time_index(i_time));
        %     topoplot(exEEG(time_samples,:), channellocations, 'MapLimits', [min(color_ticks) max(color_ticks)]);
        %     title(['Negative: ' num2str(time_index(i_time), '%.2f') ' s']);
        %     colorbar('Ticks',min(color_ticks):1:max(color_ticks));
        %     r = r + 1;
        % end
        % exEEG = mean(eegEpochs_temp.data(:,:,logical_index_positive), 3);
        % for i_time = 1:length(time_index)
        %     subplot(3,length(time_index),r)
        %     time_samples = dsearchn(params.epochTime', time_index(i_time));
        %     topoplot(exEEG(time_samples,:), channellocations, 'MapLimits', [min(color_ticks) max(color_ticks)]);
        %     title(['Positive: ' num2str(time_index(i_time), '%.2f') ' s']);
        %     r = r + 1;
        %     colorbar('Ticks',min(color_ticks):1:max(color_ticks));
        % end
        % exEEG = mean(eegEpochs_temp.data(:,:,logical_index_neutral), 3);
        % for i_time = 1:length(time_index)
        %     subplot(3,length(time_index),r)
        %     time_samples = dsearchn(params.epochTime', time_index(i_time));
        %     topoplot(exEEG(time_samples,:), channellocations, 'MapLimits', [min(color_ticks) max(color_ticks)]);
        %     title(['Neutral: ' num2str(time_index(i_time), '%.2f') ' s']);
        %     r = r + 1;
        %     colorbar('Ticks',min(color_ticks):1:max(color_ticks));
        % end
    
    % Subject-wise struct
    eegEpochs.subjectID(subID) = eegEpochs_pre;
    % Reduce data time window to 0-1s
    range_data = eegEpochs_pre.data(signal_range, :,:);
    eegEpochs.subjectID(subID).data = range_data;
end

% for subID = 1:length(subjectID)
%     eegEpochs.subjectID(subID).run_loc = run_locs{subID};
% end

%% Time frequency analysis
baseline_subtract = false;
baseline_normalise = true;
error_correct_subtraction = true;
chan = 46;
compute_theta_peak_v3_1Dcursor(eegEpochs_pre, params, baseline_subtract, baseline_normalise, error_correct_subtraction, chan);

%% ------------------------
% Organize EEG Data Per Subject
%  ------------------------

numSubjects = length(eegEpochs.subjectID); % Total subjects
allData = cell(numSubjects, 1);
allLabels = cell(numSubjects, 1);

% Put into cell arrays
for i = 1:numSubjects
    no_pos = eegEpochs.subjectID(i).data(:,:,eegEpochs.subjectID(i).label ~= 1);
    allData{i} = no_pos;
    no_pos = eegEpochs.subjectID(i).label(eegEpochs.subjectID(i).label ~= 1);
    allLabels{i} = no_pos;
end
%find_channels = {channellocations.labels}';

%% ------------------------
%  Run-Wise Cross-Validation
%  ------------------------
close all;
tuning = true; balance = false;
disp([newline 'Running fresh validation...'])
results = runWise_cv(subjectID, allData, allLabels, triggers, tuning, balance);

%% ------------------------
%  LOSO Cross-Validation
%  ------------------------
[accuracies, predictions, confusion_vals] = loso_cv(numSubjects, allData, allLabels, logspace(-6,-2,5));

%% CODE BELOW DOES NOT 

%% Build and train model for classification (single subject)

% Remove positive valence class:
data_noPos = eegEpochs_pre.data(:,:, eegEpochs_pre.label ~= 1);
label_noPos = eegEpochs_pre.label(eegEpochs_pre.label ~= 1);

% Split data between train/test
cv = cvpartition(size(data_noPos, 3), 'KFold', 5);

% Train decoder once on the full training data
trainData = data_noPos(:,:,cv.training(1));
trainLabels = label_noPos(cv.training(1));

[decoder, U1] = computeDecoder(trainData, trainLabels, 0.1);
%% Testing model
close all;
disp('Evaluate decoder on test data...');

FCB_accuracies = zeros(cv.NumTestSets, 1);
F1 = zeros(cv.NumTestSets, 1);
AUCs = zeros(numSubjects, 1);

% Test the model across folds
% figure;
for i = 1:cv.NumTestSets
    testIdx = test(cv, i);

    testData = data_noPos(:, :, testIdx);
    testLabels = label_noPos(testIdx);
    
    predictions = zeros(numel(testLabels),1);
    
    for j = 1:size(testData, 3)
        proj_test = FCB_projections(testData(:,:,testLabels == 2), testData(:,:,testLabels == 3), U1); % Apply the trained filter
        proj_test_flat = proj_test(:)';
        proj_test_sorted = sort(proj_test_flat, 'descend');
        num_selected_features = min(6, length(proj_test_sorted));
        selected_test_features = proj_test_sorted(1:num_selected_features);

        test_tbl = array2table(selected_test_features, 'VariableNames', decoder.PredictorNames);
        predictions(j) = predict(decoder, test_tbl);
    end

    FCB_accuracies(i) = sum(predictions == testLabels) / numel(testLabels);
    [~, ~, ~, AUC] = perfcurve(testLabels, predictions, 3); % Use 3 as positive class
    AUCs(i) = AUC;

    precision = sum(predictions & testLabels) / sum(predictions);
    recall = sum(predictions & testLabels) / sum(testLabels);
    F1(i) = 2 * (precision * recall) / (precision + recall);

end

disp('Visualizing performance...');
figure
plot(FCB_accuracies)
title(['Classification Accuracy Across Folds: Subject ' subjectID(6)'])
subtitle(['Features: Fisher criterion beamformer, ' 'Feature Count = ' num2str(num_selected_features) ...
    ', Regularization parameter = 0.1' newline ...
    'Classification: shrinkage-regularized LDA' ', 5-fold test split'])
xlim([0 cv.NumTestSets]); ylim([0 1])

%%

% Should I balance negative and neutral by removing neutral trials?
% Should my training and testing include data from all subjects or 5/6 and test on 1/6?


% - use weights and bias to make sigmoidal model and set threshold to .5
% - test out CCA, then power (compute_psd.m of type p_welch)
