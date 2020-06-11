function semiStruct = semiparametricfunctionReach_Force(joint,trainingInputs,trainingOutputs,testInputs,testOutputs,B1,b1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEMIPARAMETRICFUNCTION.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Descritpion: This function takes a training and test data set and finds
% the best semiparametric model to predict outputs.  It computes the outputs 
% and stores all the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15 February 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated: 
%       Derek Wolf 8 March 2019: Forces instead of torque
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPM Fit with simple model mean func and accounting for covariance of the linear
% parameters in the covariance function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate basis functions at training data
if strcmp(joint,'X')
    PHI1training = compute_phi1(trainingInputs)';
    PHI1test = compute_phi1(testInputs)';
    covSIMPLE1 = {@cov1};
    meanfunc1 = @mean1;
end
if strcmp(joint,'Y')
    PHI1training = compute_phi2(trainingInputs)';
    PHI1test = compute_phi2(testInputs)';
    covSIMPLE1 = {@cov2};
    meanfunc1 = @mean2;
end
if strcmp(joint,'Z')
    PHI1training = compute_phi3(trainingInputs)';
    PHI1test = compute_phi3(testInputs)';
    covSIMPLE1 = {@cov3};
    meanfunc1 = @mean3;
end

%B1 = diag([1000 1000]);
covGP = {@covSEardNoDer};
%covGP = {@covSEiso};
likfunc = @likGauss; 

ell = ones(size(trainingInputs,2),1)*log(0.2);
%ell = log(10);
sf = log(10);%log(0.2);

hyp1.cov = [ell;sf;B1(1,1);B1(1,2);B1(2,2)];
hyp1.lik = log(10);
hyp1.mean = b1;

covfunc1 = {'covSum',{covGP,covSIMPLE1}};

% optimize the hyperparameters and then make predictions with the optimized
% parameters
[hypSIMPLET1 fX iter] = minimize(hyp1, @gp, -100, @infExact, meanfunc1, covfunc1, likfunc, trainingInputs, trainingOutputs(:,1));
[mSIMPLET1 s2SIMPLET1] = gp(hypSIMPLET1, @infExact, meanfunc1, covfunc1, likfunc, trainingInputs, trainingOutputs(:,1), testInputs);
[mSIMPLET1train s2SIMPLET1train] = gp(hypSIMPLET1, @infExact, meanfunc1, covfunc1, likfunc, trainingInputs, trainingOutputs(:,1), trainingInputs);

% evaulate training data covariance matrix
Ky1 = covSEard(hypSIMPLET1.cov(1:size(trainingInputs,2)+1),trainingInputs)+covNoise(hypSIMPLET1.lik,trainingInputs);
K = covSEard(hypSIMPLET1.cov(1:size(trainingInputs,2)+1),trainingInputs);
Kbig = K + PHI1training*B1*PHI1training';
sn2 = exp(2*hypSIMPLET1.lik);
Kstar = covSEard(hypSIMPLET1.cov(1:size(trainingInputs,2)+1),trainingInputs,testInputs);
KstarBig = Kstar + PHI1training*B1*PHI1test';

% posterior parameters and covariance
param1Simple = inv(inv(B1)+PHI1training'*inv(Ky1)*PHI1training)*(PHI1training'*inv(Ky1)*trainingOutputs(:,1)+inv(B1)*b1);
cov1Simple = inv(inv(B1)+PHI1training'*inv(Ky1)*PHI1training);

% training error and RMS
trainingError1 = trainingOutputs(:,1)-mSIMPLET1train;
trainingRMS1 = sqrt(mean(trainingError1.^2));

% check to see if I am doing this right
R = PHI1test' - PHI1training'*inv(Ky1)*Kstar;
f = Kstar'*inv(Ky1)*trainingOutputs;
linearModel = PHI1test*param1Simple;
GPpart = Kstar'*inv(Ky1)*(trainingOutputs-PHI1training*param1Simple);
mSIMPLET2 = f + R'*param1Simple;
mSIMPLET3 = PHI1test*b1 + KstarBig'*inv(Kbig+covNoise(hypSIMPLET1.lik,trainingInputs))*(trainingOutputs-PHI1training*b1);

%%
% RMS error for RBD GPM
gpSIMPLErmsT1 = sqrt(mean((testOutputs(:,1)-mSIMPLET1).^2));
gpSIMPLErmsT2 = sqrt(mean((testOutputs(:,1)-mSIMPLET2).^2));
gpSIMPLErmsT3 = sqrt(mean((testOutputs(:,1)-mSIMPLET3).^2));

gpSIMPLEErrors = testOutputs - [mSIMPLET1];
gpSIMPLErSqr = 1 - sum(gpSIMPLEErrors.^2)./var(testOutputs)/length(testOutputs);

% compute the integrated variance in the predictions for all test samples
intVar1 = sum(s2SIMPLET1);

%close all
h5 = figure(5);
z = linspace(1,length(mSIMPLET1),length(mSIMPLET1));
hold on
hv = fill([z'; flipdim(z',1)],[mSIMPLET1+2*sqrt(s2SIMPLET1);flipdim(mSIMPLET1-2*sqrt(s2SIMPLET1),1)],[7 7 7]/8);
htSIMPLET1 =plot(testOutputs(:,1),'b','LineWidth',2);
hmSIMPLET1 = plot(mSIMPLET1,'r--');
hmSIMPLET4 = plot(GPpart,'c');
hmSIMPLET5 = plot(linearModel,'k');
legend([htSIMPLET1 hmSIMPLET1 hmSIMPLET4 hmSIMPLET5],'actual','predicted','GP','LM')
title([joint,' torque'])
ylabel('torque (N-m)')
saveas(h5,['./Modelfigures/',joint,'Semiparametric.jpg'])
saveas(h5,['./Modelfigures/',joint,'Semiparametric.fig'])
saveas(h5,['./Modelfigures/',joint,'Semiparametric.eps'],'epsc')
hold off

% build the structure of variables to save
semiStruct.trainingError        = trainingError1;
semiStruct.trainingRMS          = trainingRMS1;
semiStruct.trainingPredictions  = mSIMPLET1train;
semiStruct.trainingInputs       = trainingInputs;
semiStruct.trainingOutputs      = trainingOutputs;
semiStruct.testRMS              = gpSIMPLErmsT1;
semiStruct.testRsqr             = gpSIMPLErSqr;
semiStruct.testInputs           = testInputs;
semiStruct.testOutputs          = testOutputs;
semiStruct.testErrors           = gpSIMPLEErrors;
semiStruct.testVariance         = s2SIMPLET1;
semiStruct.testPredictions      = mSIMPLET1;
semiStruct.priorParameters      = b1;
semiStruct.priorParameterCov    = B1;
semiStruct.posteriorParameters  = param1Simple;
semiStruct.posteriorParameterCov= cov1Simple;
semiStruct.gpHyperparameters    = hypSIMPLET1;

