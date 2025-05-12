function [predictions, posteriors] = classifyTestData(testData, decoder, U1, num_channels)
    num_trials = size(testData, 3);
    predictions = zeros(num_trials, 1);
    posteriors = zeros(num_trials, 2);

    testData = permute(testData, [2 1 3]);

    for j = 1:num_trials
        %proj_test = U1' * squeeze(testData(:,:,j)); % Project test trial
        proj_test = squeeze(testData(:,:,j))' * U1; % Project test trial

        reduced_test = proj_test(:,1:num_channels);

        % Visualize discrimancy (figure/output 1-3 are negative, 4-5 are neutral
        % if j < 4
        %     figure
        %     plot(proj_test(:,1)); hold on; plot(proj_test(:,2)); hold on; plot(proj_test(:,3)); % Compare across classes
        %     disp(mean(mean(reduced_test,1)))
        % elseif (j > num_trials - 6 && j < num_trials - 3)
        %     figure
        %     plot(proj_test(:,1)); hold on; plot(proj_test(:,2)); hold on; plot(proj_test(:,3)); % Compare across classes
        %     disp(mean(mean(reduced_test,1)))
        % end

        downsample_rate = 50;
        reduced_test = reduced_test(1:downsample_rate:size(reduced_test,1),:);
        
        selected_test_features = reduced_test(:)'; % proj_test_flat = proj_test(:)';
        % proj_test_sorted = sort(proj_test_flat, 'descend');
        % num_selected_features = min(10, length(proj_test_sorted));
        % selected_test_features = proj_test_sorted(1:num_selected_features);

        test_tbl = array2table(selected_test_features, 'VariableNames', decoder.PredictorNames);
        [predictions(j), posteriors(j,:)] = predict(decoder, test_tbl);
    end
end