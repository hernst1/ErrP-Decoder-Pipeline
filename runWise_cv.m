function results = runWise_cv(subjectID, allData, allLabels, triggers, tuning, balance)

results = struct;
numSubjects = size(allData,1);

% tuning = 1;
% balance = 1;
num_channels = 3;

for sub = 1:numSubjects

    disp(['Subject ' num2str(subjectID(sub))])

    data = allData{sub};
    labels = allLabels{sub};

    % Extract labels from trigger_order
    idx = triggers.subjectID(sub).ordered(:,3) ~= 1;
    ordered_mat = triggers.subjectID(sub).ordered(idx,:);
    ordered_vec = ordered_mat(:, 2);
    %ordered = ordered(triggers.subjectID(sub).ordered(:, 3) ~= 1);
    idx = triggers.subjectID(sub).unordered(:,2) ~= 1;
    unordered_mat = triggers.subjectID(sub).unordered(idx,:);
    unordered_vec = unordered_mat(:, 1);
    %unordered = unordered(triggers.subjectID(sub).unordered(:, 2) ~= 1);
    [~, sort_idx] = ismember(unordered_vec, ordered_vec);
    sorted_triggers = ordered_mat(sort_idx,:);
    runs = sorted_triggers(:,1);

    numRuns = length(unique(runs));
    for testRun = 1:numRuns
        disp(['Leaving out Run ' num2str(testRun) ' from training set'])
        trainIdx = runs ~= testRun;

        rng(2);
        % Even out classes by sampling from testData
        testData = data(:,:,~trainIdx);
        testLabels = labels(~trainIdx);
        if balance
            min_class = min(sum(testLabels == 2), sum(testLabels == 3));
            disp(['Testing on ' num2str(min_class) ' from each class...'])
            sampled_testIdx = [randsample(find(testLabels == 2), min_class); randsample(find(testLabels == 3), min_class)];
            sampled_testData = testData(:,:,sampled_testIdx);
            sampled_testLabels = testLabels(sampled_testIdx);
            % disp(['Testing data has ' num2str(sum(sampled_testLabels == 2)) ' samples from the negative class' newline ...
            %     ' and ' num2str(sum(sampled_testLabels == 3)) ' samples from the neutral class.'])

            testData = sampled_testData;
            testLabels = sampled_testLabels;
        end

        trainData = data(:,:,trainIdx);
        trainLabels = labels(trainIdx);

        if balance
            min_class = min(sum(trainLabels == 2), sum(trainLabels == 3));
            disp(['Training on ' num2str(min_class) ' from each class...'])
            
            sampled_trainIdx = [randsample(find(trainLabels == 2), min_class); randsample(find(trainLabels == 3), min_class)];
            sampled_trainData = trainData(:,:,sampled_trainIdx);
            sampled_trainLabels = trainLabels(sampled_trainIdx);
            % disp(['Training data has ' num2str(sum(sampled_trainLabels == 2)) ' samples from the negative class' newline ...
            %     ' and ' num2str(sum(sampled_trainLabels == 3)) ' samples from the neutral class.'])

            trainData = sampled_trainData;
            trainLabels = sampled_trainLabels;
        end

        % ------------------------
        %  Hyperparameter (Gamma, γ) Tuning
        %  ------------------------
        best_gamma = 10e-3; % Default regularization value
        if tuning
            best_acc = 0;
            gamma_values = logspace(-7,-2,6);
            for gamma = gamma_values
                %disp(['Testing γ = ', num2str(gamma)]);
                [decoder, U1] = computeDecoder(trainData, trainLabels, 'CCA', gamma);
    
                val_predictions = classifyTestData(testData, decoder, U1, num_channels);
                valLabels = testLabels;
                val_acc = sum(val_predictions == valLabels) / numel(valLabels);
                
                if val_acc >= best_acc
                    best_acc = val_acc;
                    best_gamma = gamma;
                end
            end
        end
        disp(['Using γ = ', num2str(best_gamma), ' for Run ', num2str(testRun)]);
        [decoder, U1] = computeDecoder(trainData, trainLabels, 'CCA', best_gamma);
        
        % plot posteriors by class or use perfcurve to get auc
        [predictions, posteriors] = classifyTestData(testData, decoder, U1, num_channels);
        accuracies = sum(predictions == testLabels) / numel(testLabels);
        [X, Y, T, AUC] = perfcurve(testLabels, posteriors(:,decoder.ClassNames == 2), 2); % Use 2 as error (negative) class
        % figure
        % plot(X,Y)

        C = confusionmat(testLabels,predictions);

        tpr = C(1,1)/sum(C(1,:));
        tnr = C(2,2)/sum(C(2,:));
        confusion_vals = [tpr 1-tpr; 1-tnr tnr];
        makeConfusion(confusion_vals, decoder)

        results.subject(sub).run(testRun).gamma = best_gamma;
        results.subject(sub).run(testRun).confusion = confusion_vals;
        results.subject(sub).run(testRun).accuracies = accuracies;
        results.subject(sub).run(testRun).posteriors = posteriors;
        results.subject(sub).run(testRun).AUC = AUC;
        %disp(selected_channels)
    end
end

end