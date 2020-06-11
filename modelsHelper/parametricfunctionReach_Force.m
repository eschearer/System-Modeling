function paraStruct = parametricfunctionReach_Force(joint,trainingInputs,trainingOutputs,testInputs,testOutputs,Bprior,bprior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRICFUNCTION.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Descritpion: This function takes a training and test data set and finds
% the best parametric model to predict outputs.  It computes the outputs 
% and stores all the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15 February 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametric model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate basis functions at training and test data
if strcmp(joint,'X')
    PHI1training    = compute_phi1(trainingInputs)';
    PHI1test        = compute_phi1(testInputs)';
end
if strcmp(joint,'Y')
    PHI1training    = compute_phi2(trainingInputs)';
    PHI1test        = compute_phi2(testInputs)';
end
if strcmp(joint,'Z')
    PHI1training    = compute_phi3(trainingInputs)';
    PHI1test        = compute_phi3(testInputs)';
end

% compute the model parameters from the training data
beta1 = pinv(PHI1training)*trainingOutputs(:,1);
% predictions of training outputs
simpleTrainingPred1 = PHI1training*beta1;
% training error
trainingError1 = simpleTrainingPred1 - trainingOutputs(:,1);
% compute the estimated noise from the residuals of the training data
sqResid1 = mean((simpleTrainingPred1 - trainingOutputs(:,1)).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now do it again using Bayesian regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B1 = inv(PHI1training'*PHI1training/sqResid1+inv(Bprior));
beta1 = B1*(inv(Bprior)*bprior'+PHI1training'*trainingOutputs(:,1)/sqResid1);
% predictions of training outputs
simpleTrainingPred1 = PHI1training*beta1;
% training error
trainingError1 = simpleTrainingPred1 - trainingOutputs(:,1);
% compute the estimated noise from the residuals of the training data
sqResid1 = mean((simpleTrainingPred1 - trainingOutputs(:,1)).^2);
% training RMS
trainingRMS1 = sqrt(sqResid1);

%%
% predict the test outputs

simplePred1 = PHI1test*beta1;
% compute the variance of the test predictions
simpleVar1 = zeros(length(simplePred1),1);
for i = 1:length(simplePred1)
    simpleVar1(i) = sqResid1 + (PHI1test(i,:)*B1*PHI1test(i,:)');
end

% compute the information matrix of the training data, its determinant, and
% the integrated variance in the predictions for all test samples
info1 = PHI1training'*PHI1training/sqResid1+inv(Bprior);
det1 = det(info1);
intVar1 = sum(simpleVar1);

% RMS error for 
simplermsT1 = sqrt(mean((testOutputs(:,1)-simplePred1).^2));

simpleErrors = testOutputs - [simplePred1];
simplerSqr = 1 - sum(simpleErrors.^2)./var(testOutputs)/length(testOutputs);

close all
h9 = figure(9);
z = linspace(1,length(simplePred1),length(simplePred1));
hold on
hv = fill([z'; flipdim(z',1)],[simplePred1+2*sqrt(simpleVar1);flipdim(simplePred1-2*sqrt(simpleVar1),1)],[7 7 7]/8);
htSIMPLET1 =plot(testOutputs(:,1),'b','LineWidth',2);
hsimplePred1 = plot(simplePred1,'r--');
legend([htSIMPLET1 hsimplePred1],'actual','predicted')
title([joint,' torque'])
ylabel('torque (N-m)')
saveas(h9,['./Modelfigures/',joint,'Parametric.jpg'])
saveas(h9,['./Modelfigures/',joint,'Parametric.fig'])
saveas(h9,['./Modelfigures/',joint,'Parametric.eps'],'epsc')
hold off

% build the structure of variables to save
paraStruct.trainingError        = trainingError1;
paraStruct.trainingRMS          = trainingRMS1;
paraStruct.trainingPredictions  = simpleTrainingPred1;
paraStruct.trainingInputs       = trainingInputs;
paraStruct.trainingOutputs      = trainingOutputs;
paraStruct.infomationMatrix     = info1;
paraStruct.testRMS              = simplermsT1;
paraStruct.testRsqr             = simplerSqr;
paraStruct.testInputs           = testInputs;
paraStruct.testOutputs          = testOutputs;
paraStruct.testErrors           = simpleErrors;
paraStruct.testVariance         = simpleVar1;
paraStruct.testPredictions      = simplePred1;
paraStruct.parameters           = beta1;
paraStruct.parameterCov         = B1;
