%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: computeDecoder.m
%
% Function Description: computeDecoder.m extracts fisher criterion
% beamformer features for each trial in the training dataset. The spatial 
% filter is projected on the dataset
% 
% Function Usage Instructions: The function is called as follows:
% >> [decoder, V1] = computeDecoder(trainEpochs, trainLabels);
%
% Inputs:
%   trainEpochs: training data with shape samples x channels x trials
%   trainLabels: training labels with shape labels x 1
%
% Outputs:
%   decoder: classification discriminant class decoder
%   V1: eigenvalues from FCB spatial filter
%
% Example usage:
% >> [decoder, V1] = computeDecoder(trainData, testData);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [decoder, U1] = computeDecoder(trainEpochs, trainLabels, feature_type, gamma)
%%%%%%%%%%%%%%%%%%
%% Feature Extraction using Fisher Criterion Beamformer (FCB) %%
%%%%%%%%%%%%%%%%%%%% 
%disp('Extracting features using FCB...');

num_samples = size(trainEpochs, 1);
num_channels = size(trainEpochs, 2);
num_trials = size(trainEpochs, 3);

% Compute class-wise mean and covariance
X_neutral = trainEpochs(:, :, trainLabels == 3);
X_negative = trainEpochs(:, :, trainLabels == 2);

switch(feature_type)
    case 'FCB'
        [U1, ~] = FCB_spatial_filters(X_negative, X_neutral, gamma);
        [X_negative_f, X_neutral_f] = FCB_projections(X_negative, X_neutral, U1);
    
        select_channels = 3;
        X_negative_f = X_negative_f(1:select_channels,:,:);
        X_neutral_f = X_neutral_f(1:select_channels,:,:);
        %selected_channels = sorted(1:select_channels);
        
        downsample_rate = 50;
        X_negative_f_downsampled = X_negative_f(:,1:downsample_rate:length(X_negative_f),:);
        X_neutral_f_downsampled = X_neutral_f(:,1:downsample_rate:length(X_neutral_f),:);
        
        num_features = size(X_negative_f_downsampled,2) * select_channels;
        features = zeros(num_trials, num_features);
        
        neut_count = 0;
        neg_count = 0;
        
        for i = 1:num_trials
            if trainLabels(i) == 3
                neut_count = neut_count + 1;
                proj = squeeze(X_neutral_f_downsampled(:,:,neut_count));
            elseif trainLabels(i) == 2
                neg_count = neg_count + 1;
                proj = squeeze(X_negative_f_downsampled(:,:,neg_count));
            else
                continue;
            end
            features(i, :) = proj(:)';
        end
        selected_features_tbl = array2table(features);

    case 'CCA'
        U1 = get_cca_spatialfilter(trainEpochs, trainLabels);
        
        select_channels = 3;
        downsample_rate = 50;

        for i = 1:num_trials
            proj = squeeze(trainEpochs(:,:,i)) * U1;
            
            X = proj(:,1:select_channels);
            X_downsampled = X(1:downsample_rate:length(X),:);
        
            features(i, :) = X_downsampled(:)';
        end
        selected_features_tbl = array2table(features);
end

%%%%%%%%%%%%%%%%%%
%% Train Model
%%%%%%%%%%%%%%%%%%%% 
%disp('Training classifier using sLDA...');
% LDA
% Compute class weights
classes = unique(trainLabels);
counts = histc(trainLabels, classes);
weights = 1 ./ counts;
sample_weights = weights(arrayfun(@(x) find(classes == x), trainLabels));

%prior = struct('ClassNames', unique(trainLabels), 'ClassProbs', [0.75, 0.25]);  % or equal weights if balanced


decoder = fitcdiscr(selected_features_tbl, trainLabels, 'DiscrimType', 'linear', 'Gamma', gamma, 'Weights', sample_weights, 'Prior','uniform');%, 'Prior', prior);%
% SVM
%decoder = fitcsvm(selected_features_tbl, trainLabels, 'Standardize', false, 'KernelFunction', 'gauss', 'KernelScale', 'Auto');%, 'Prior', prior);

end