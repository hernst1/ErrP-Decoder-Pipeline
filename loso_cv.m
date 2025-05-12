%% ------------------------
function [accuracies, predictions, confusion_vals] = loso_cv(numSubjects, allData, allLabels, gamma_values)
%  ------------------------

accuracies = zeros(numSubjects, 1);
for testSubj = 1:numSubjects
    %testSubj = 1;
    disp(['Cross-validation Fold ', num2str(testSubj), ' / ', num2str(numSubjects)]);

    % Split Training & Testing Data by Subject
    trainData = [];
    trainLabels = [];
    for i = 1:numSubjects
        if i ~= testSubj  % Exclude test subject
            trainData = cat(3, trainData, allData{i});
            trainLabels = [trainLabels; allLabels{i}];
        end
    end
    testData = allData{testSubj};
    testLabels = allLabels{testSubj};

    % ----------------------
    %  Hyperparameter Tuning (Optimize γ for FCB Regularization)
    % ----------------------
    best_gamma = 0.1; % Default regularization value
    best_acc = 0;

    %gamma_values = [0.00001];
    %gamma_values = logspace(-6, -3, 4); % Test γ from 10^-4 to 1

    for gamma = gamma_values
        disp(['Testing γ = ', num2str(gamma)]);
        [decoder, U1] = computeDecoder(trainData, trainLabels, 'FCB', gamma);
        
        % Validation split within training data
        cv = cvpartition(length(trainLabels), 'HoldOut', 0.2);
        valData = trainData(:,:,cv.test);
        valLabels = trainLabels(cv.test);
        
        val_predictions = classifyTestData(valData, decoder, U1);
        val_acc = sum(val_predictions == valLabels) / numel(valLabels);
        
        if val_acc > best_acc
            best_acc = val_acc;
            best_gamma = gamma;
        end
    end
    
    % ----------------------
    %  Train Final Model Using Best γ (too expensive to test all)
    % ----------------------
    disp(['Using γ = ', num2str(best_gamma), ' for Subject ', num2str(testSubj)]);
    [decoder, U1] = computeDecoder(trainData, trainLabels, 'FCB', best_gamma);
    
    % ----------------------
    %  Test Model on Left-Out Subject
    % ----------------------
    predictions = classifyTestData(testData, decoder, U1);
    accuracies(testSubj) = sum(predictions == testLabels) / numel(testLabels);
    
    % Compute Confusion
    C = confusionmat(testLabels,predictions);
    tpr = C(1,1)/sum(C(1,:));
    tnr = C(2,2)/sum(C(2,:));
    confusion_vals = [tpr 1-tpr; 1-tnr tnr];
    makeConfusion(confusion_vals, decoder)
end

% ------------------------
%  Report Final Performance Metrics
%  ------------------------
disp(['Mean Accuracy: ', num2str(mean(accuracies) * 100), ' ± ', num2str(std(accuracies) * 100)]);



