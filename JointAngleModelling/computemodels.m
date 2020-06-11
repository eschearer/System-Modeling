%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTEMODELS.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Descritpion: This script loads the 'parseddata.mat' file and does GPR for
% each muscle where the inputs are joint angles and the outputs are joint
% torques.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15 February 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

load('parseddata')

% prior parameters and covariance
params1 = zeros(1,2);
params2 = zeros(1,2);
params3 = zeros(1,2);
params4 = zeros(1,2);
covariance1 = diag(10000*ones(length(params1),1));
covariance2 = diag(10000*ones(length(params2),1));
covariance3 = diag(10000*ones(length(params3),1));
covariance4 = diag(10000*ones(length(params4),1));

muscles = 10;
targets = 27;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build data structure to save model data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modeldata.muscle(1).label  = 'Radial Nerve - Triceps';
modeldata.muscle(2).label  = 'Axillary Nerve - Deltoids';
modeldata.muscle(3).label  = 'Thoracodorsal Nerve - Latissimus Dorsi';
modeldata.muscle(4).label  = 'Long Thoracic - Serratus Anterior';
modeldata.muscle(5).label  = 'Musculocutaneous Nerve - Biceps/Bracialis';
modeldata.muscle(6).label  = 'Suprascapular Nerve - Supraspinatus/Infraspinatus';
modeldata.muscle(7).label  = 'Rhomboids';
modeldata.muscle(8).label  = 'Lower Pectoralis';
modeldata.muscle(9).label  = 'Upper Pectoralis';
modeldata.muscle(10).label = 'Passive';

for i = 1:muscles
    i
    angleData   = [];
    torqueData  = [];
    for j = 1:targets
        if size(stimdata.muscle(i).target(j).allJointTorques,1) > 0
            torqueData  = [torqueData;stimdata.muscle(i).target(j).allJointTorques];
            angleData   = [angleData;ones(size(stimdata.muscle(i).target(j).allJointTorques,1),1)*[stimdata.muscle(i).target(j).meanFilteredJointAngles(1:4) stimdata.muscle(i).target(j).meanFilteredJointAngles(6)]];
        end
    end
    modeldata.muscle(i).elevationplane.parametric           = parametricfunctionReach('ElevationPlane',   angleData,torqueData(:,1),angleData,torqueData(:,1),covariance1,params1);
    modeldata.muscle(i).elevationplane.semiparametric       = semiparametricfunctionReach('ElevationPlane',   angleData(1:2:length(angleData),:),torqueData(1:2:length(angleData),1),angleData,torqueData(:,1),modeldata.muscle(i).elevationplane.parametric.parameterCov,   modeldata.muscle(i).elevationplane.parametric.parameters);
    %
    display('a');
    modeldata.muscle(i).shoulderelevation.parametric        = parametricfunctionReach('ShoulderElevation',angleData,torqueData(:,2),angleData,torqueData(:,2),covariance1,params1);
    modeldata.muscle(i).shoulderelevation.semiparametric    = semiparametricfunctionReach('ShoulderElevation',angleData(1:2:length(angleData),:),torqueData(1:2:length(angleData),2),angleData,torqueData(:,2),modeldata.muscle(i).shoulderelevation.parametric.parameterCov,modeldata.muscle(i).shoulderelevation.parametric.parameters);
    %
    display('b');
    modeldata.muscle(i).shoulderrotation.parametric         = parametricfunctionReach('ShoulderRotation', angleData,torqueData(:,3),angleData,torqueData(:,3),covariance1,params1);
    modeldata.muscle(i).shoulderrotation.semiparametric     = semiparametricfunctionReach('ShoulderRotation', angleData(1:2:length(angleData),:),torqueData(1:2:length(angleData),3),angleData,torqueData(:,3),modeldata.muscle(i).shoulderrotation.parametric.parameterCov, modeldata.muscle(i).shoulderrotation.parametric.parameters);
    %
    display('c');
    modeldata.muscle(i).elbowflexion.parametric             = parametricfunctionReach('ElbowFlexion',     angleData,torqueData(:,4),angleData,torqueData(:,4),covariance1,params1);
    modeldata.muscle(i).elbowflexion.semiparametric         = semiparametricfunctionReach('ElbowFlexion',     angleData(1:2:length(angleData),:),torqueData(1:2:length(angleData),4),angleData,torqueData(:,4),modeldata.muscle(i).elbowflexion.parametric.parameterCov,     modeldata.muscle(i).elbowflexion.parametric.parameters);
    %
    display('d');
    modeldata.muscle(i).elbowpronation.parametric           = parametricfunctionReach('ElbowPronation',   angleData,torqueData(:,5),angleData,torqueData(:,5),covariance1,params1);
    modeldata.muscle(i).elbowpronation.semiparametric       = semiparametricfunctionReach('ElbowPronation',   angleData(1:2:length(angleData),:),torqueData(1:2:length(angleData),5),angleData,torqueData(:,5),modeldata.muscle(i).elbowpronation.parametric.parameterCov,   modeldata.muscle(i).elbowpronation.parametric.parameters);
    %
end
        
save('./modelfiles/models30May2017.mat','modeldata')

