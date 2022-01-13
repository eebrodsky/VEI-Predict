function [trainedClassifier, validationAccuracy] = trainClassifier_Generic(trainingData, Predictors)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% Returns a trained classifier and its accuracy. This code recreates the
% classification model trained in Classification Learner app. Use the
% generated code to automate training the same model with new data, or to
% learn how to programmatically train models.
%
%  Input:
%      trainingData: A table containing the same predictor and response
%       columns as those imported into the app.
%
%  Output:
%      trainedClassifier: A struct containing the trained classifier. The
%       struct contains various fields with information about the trained
%       classifier.
%
%      trainedClassifier.predictFcn: A function to make predictions on new
%       data.
%
%      validationAccuracy: A double containing the accuracy in percent. In
%       the app, the History list displays this overall accuracy score for
%       each model.
%
% Use the code to train the model with new data. To retrain your
% classifier, call the function from the command line with your original
% data or new data as the input argument trainingData.
%
% For example, to retrain a classifier trained with the original data set
% T, enter:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% To make predictions with the returned 'trainedClassifier' on new data T2,
% use
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 must be a table containing at least the same predictor columns as used
% during training. For details, enter:
%   trainedClassifier.HowToPredict
% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictorNames = Predictors; %{'medianVEI', 'morphology', 'RockType', 'TectonicSetting', 'meanVEI', 'minVEI', 'maxVEI', 'modeVEI', 'eruptiontimediff', 'eruptionduration'};
predictors = inputTable(:, predictorNames);
response = inputTable.currentVEI;
%isCategoricalPredictor = [false, true, true, true, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
% template = templateSVM(...
%     'KernelFunction', 'gaussian', ...
%     'PolynomialOrder', [], ...
%     'KernelScale', 3.2, ... 
%     'BoxConstraint', 1, ...
%     'Standardize', true);
% classificationSVM = fitcecoc(...
%     predictors, ...
%     response, ...
%     'Learners', template, ...
%     'Coding', 'onevsone');

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateTree('MaxNumSplits',5, ...
    'NumVariablesToSample','all'); 
n=height(inputTable);
kf=5; % kfold partition

classificationEnsemble = fitcensemble(...
    predictors, ...
    response, ...
    'Learners',template, ...
    'OptimizeHyperparameters', {'NumLearningCycles','Method','LearnRate'}, ....
    'HyperparameterOptimizationOptions', ...
    struct('Verbose',1,'ShowPlots',false,'KFold',kf,'UseParallel',1));
 

% 'OptimizeHyperparameters', {'NumLearningCycles','MaxNumSplits','Method','LearnRate'},
% Create the result struct with predict function
% predictorExtractionFcn = @(t) t(:, predictorNames);
% svmPredictFcn = @(x) predict(classificationSVM, x);
% trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));
predictorExtractionFcn = @(t) t(:, predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
%trainedClassifier.ClassificationSVM = classificationSVM;
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
% 
 rng('default'); % for reproducibility
 NValidationCalculations=25;
% 
parfor i=1:NValidationCalculations,
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', kf);

% Compute validation predictions
%[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
va(i) = 1 - kfoldLoss(partitionedModel,'LossFun', 'ClassifError');
%p=round(kfoldPredict(partitionedModel));
%va(i)=mean(response==p);
end
% 
 validationAccuracy = mean(va);