%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computeModels_wrist.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function computemodels_wristPosition_global(folder)
close all


load([folder,'/parseddata'])

filename=['model_wristPosition_',folder,'.mat'];

params = [1 1 1 1]*log(10);

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
    wristData=[];
    forceData=[];
    for j=1:targets
        if size(stimdata.muscle(i).target(j).allMeanNewHmForce,1)>0
            if ~isnan(stimdata.muscle(10).target(j).meanWristPosition_global(:,1))
                forceData=[forceData;stimdata.muscle(i).target(j).allMeanNewHmForce];
                wristData=[wristData;ones(size(stimdata.muscle(i).target(j).allMeanNewHmForce,1),1)*[stimdata.muscle(i).target(j).meanWristPosition_global]];
            end
        end
    end
    testOutputs = forceData;
    testInputs = wristData;
    
    modeldata.muscle(i).Fx = TrainWristGP_position('Fx',i,wristData,forceData(:,1),params);
    modeldata.muscle(i).Fy = TrainWristGP_position('Fy',i,wristData,forceData(:,2),params);
    modeldata.muscle(i).Fz = TrainWristGP_position('Fz',i,wristData,forceData(:,3),params);
    
    
end

save(['models_global/',filename],'modeldata');
