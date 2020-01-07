%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GETACTIVE.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Descritpion: This script opens data files, computes joint angles with RTS
% smoother, finds the portions where the muscles were fully activated, and 
% sorts the data by muscle, target, and rep.  It saves a grand data file 
% called 'parseddata.mat'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 20 November 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated:
% 29 Aug 17:    Changed the time to cut off at the start of each trial from
%               0.25 seconds to 0.5 seconds.  Added means for each trial so
%               we don't capture 13 hz variation in force in our GP models.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%
% define time parameters
%%%%%%%%%%%%%%%%%%%%%%%%
stimTime = 1.25;
restTime = 0.75;
riseTime = 0.5; % changed this from 0.25 on 29 Aug 17
global dt;
dt = 1/52;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the data structure for storing things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muscles = 10;
targets = 27;
reps = 3;
for i = 1:muscles
    for j = 1:targets
        for k = 1:reps
            stimdata.muscle(i).target(j).rep(k).handMarker          = zeros(1/dt:4);
            stimdata.muscle(i).target(j).rep(k).hmForce             = zeros(1/dt:3);
            stimdata.muscle(i).target(j).rep(k).newHmForce          = zeros(1/dt:3);
            stimdata.muscle(i).target(j).rep(k).hmMeasPosition      = zeros(1/dt:3);
            stimdata.muscle(i).target(j).rep(k).stimLevel           = zeros(1/dt:1);
            stimdata.muscle(i).target(j).rep(k).stimPW              = zeros(1/dt:12);
            stimdata.muscle(i).target(j).rep(k).thoraxQuaternion    = zeros(1/dt:4);
            stimdata.muscle(i).target(j).rep(k).wristPosition       = zeros(1/dt:3);
            stimdata.muscle(i).target(j).rep(k).wristPosition_global       = zeros(1/dt:3);
            stimdata.muscle(i).target(j).rep(k).wristQuaternion     = zeros(1/dt:4);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% load the target order
%%%%%%%%%%%%%%%%%%%%%%%
load('targetOrder1.txt')
load('targetOrder2.txt')
load('targetOrder3.txt')

targetOrder = [targetOrder1 targetOrder2 targetOrder3];
repOrder = [ones(targets,1) 2*ones(targets,1) 3*ones(targets,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through all of the data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trials = targets*reps;
trialIndex = 1:trials; 
for l = 1:length(trialIndex)
    %close all;
    trial = trialIndex(l);
    load(['dataTrial',num2str(trial)])
    trial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % use RTS smoother to find joint angles 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    smoothmotion_wrist
    %pause 
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find where activation starts for each muscle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lastMuscle = 0;
    startIndex = []; % vector of the index each time stimulation begins
    for i = 1:length(time)
        if muscle(i) ~= 0 && lastMuscle == 0 % first time a muscle is activated
            startIndex = [startIndex i];
            lastMuscle = 1;
        else
            lastMuscle = muscle(i);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % record data for the last 1 second of each stim
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currentTarget   = targetOrder(trial);
    currentRep      = repOrder(trial);
    for i = 1:length(startIndex)
        
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).handMarker         = handMarker(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).hmForce            = hmForce(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).newHmForce         = newHmForce(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).hmMeasPosition     = hmMeasPosition(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).stimLevel          = stimLevel(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).stimPW             = stimPW(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).thoraxQuaternion   = thoraxQuaternion(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).wristPosition      = filteredWristPosition(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).wristPosition_global= filteredWristPosition_global(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).wristQuaternion    = forearmQuaternion_thor(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        
        % these mean values so that we are not capturing extra noise in the models
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).meanHandMarker         = nanmean(stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).handMarker); 
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).meanHmForce            = nanmean(stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).hmForce); 
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).meanNewHmForce         = nanmean(stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).newHmForce); 
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).meanHmMeasPosition     = nanmean(stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).hmMeasPosition); 
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).meanThoraxQuaternion   = nanmean(stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).thoraxQuaternion); 
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).meanWristPosition      = nanmean(stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).wristPosition); 
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).meanWristPosition_global= nanmean(stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).wristPosition_global); 
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).meanWristQuaternion    = nanmean(stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).wristQuaternion); 

        
    end
    %pause
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each muscle and target put all reps into one bin and average joint
% angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:muscles
    for j = 1:targets
        stimdata.muscle(i).target(j).allHandMarker          = [];
        stimdata.muscle(i).target(j).allHmForce             = [];
        stimdata.muscle(i).target(j).allNewHmForce          = [];
        stimdata.muscle(i).target(j).allHmMeasPosition      = [];
        stimdata.muscle(i).target(j).allStimLevel           = [];
        stimdata.muscle(i).target(j).allStimPW              = [];
        stimdata.muscle(i).target(j).allThoraxQuaternion    = [];
        stimdata.muscle(i).target(j).allWristPosition       = [];
        stimdata.muscle(i).target(j).allWristPosition_global= [];
        stimdata.muscle(i).target(j).allWristQuaternion     = [];
        
        stimdata.muscle(i).target(j).allMeanHandMarker          = [];
        stimdata.muscle(i).target(j).allMeanHmForce             = [];
        stimdata.muscle(i).target(j).allMeanNewHmForce          = [];
        stimdata.muscle(i).target(j).allMeanHmMeasPosition      = [];
        stimdata.muscle(i).target(j).allMeanThoraxQuaternion    = [];
        stimdata.muscle(i).target(j).allMeanWristPosition       = [];
        stimdata.muscle(i).target(j).allMeanWristPosition_global= [];
        stimdata.muscle(i).target(j).allMeanWristQuaternion     = [];
        stimdata.muscle(i).target(j).allMeanNewHmForce_Thor     = []; %Force in thorax Frame
        
        for k = 1:reps
            stimdata.muscle(i).target(j).allHandMarker          = [stimdata.muscle(i).target(j).allHandMarker; stimdata.muscle(i).target(j).rep(k).handMarker];
            stimdata.muscle(i).target(j).allHmForce             = [stimdata.muscle(i).target(j).allHmForce; stimdata.muscle(i).target(j).rep(k).hmForce];
            stimdata.muscle(i).target(j).allNewHmForce          = [stimdata.muscle(i).target(j).allNewHmForce; stimdata.muscle(i).target(j).rep(k).newHmForce];
            stimdata.muscle(i).target(j).allHmMeasPosition      = [stimdata.muscle(i).target(j).allHmMeasPosition; stimdata.muscle(i).target(j).rep(k).hmMeasPosition];
            stimdata.muscle(i).target(j).allStimLevel           = [stimdata.muscle(i).target(j).allStimLevel; stimdata.muscle(i).target(j).rep(k).stimLevel];
            stimdata.muscle(i).target(j).allStimPW              = [stimdata.muscle(i).target(j).allStimPW; stimdata.muscle(i).target(j).rep(k).stimPW];
            stimdata.muscle(i).target(j).allThoraxQuaternion    = [stimdata.muscle(i).target(j).allThoraxQuaternion; stimdata.muscle(i).target(j).rep(k).thoraxQuaternion];
            stimdata.muscle(i).target(j).allWristPosition       = [stimdata.muscle(i).target(j).allWristPosition; stimdata.muscle(i).target(j).rep(k).wristPosition];
            stimdata.muscle(i).target(j).allWristPosition_global= [stimdata.muscle(i).target(j).allWristPosition_global; stimdata.muscle(i).target(j).rep(k).wristPosition_global];
            stimdata.muscle(i).target(j).allWristQuaternion     = [stimdata.muscle(i).target(j).allWristQuaternion; stimdata.muscle(i).target(j).rep(k).wristQuaternion];
            
            stimdata.muscle(i).target(j).allMeanHandMarker          = [stimdata.muscle(i).target(j).allMeanHandMarker; stimdata.muscle(i).target(j).rep(k).meanHandMarker];
            stimdata.muscle(i).target(j).allMeanHmForce             = [stimdata.muscle(i).target(j).allMeanHmForce; stimdata.muscle(i).target(j).rep(k).meanHmForce];
            stimdata.muscle(i).target(j).allMeanNewHmForce          = [stimdata.muscle(i).target(j).allMeanNewHmForce; stimdata.muscle(i).target(j).rep(k).meanNewHmForce];
            stimdata.muscle(i).target(j).allMeanHmMeasPosition      = [stimdata.muscle(i).target(j).allMeanHmMeasPosition; stimdata.muscle(i).target(j).rep(k).meanHmMeasPosition];
            stimdata.muscle(i).target(j).allMeanThoraxQuaternion    = [stimdata.muscle(i).target(j).allMeanThoraxQuaternion; stimdata.muscle(i).target(j).rep(k).meanThoraxQuaternion];
            stimdata.muscle(i).target(j).allMeanWristPosition       = [stimdata.muscle(i).target(j).allMeanWristPosition; stimdata.muscle(i).target(j).rep(k).meanWristPosition];
            stimdata.muscle(i).target(j).allMeanWristPosition_global= [stimdata.muscle(i).target(j).allMeanWristPosition_global; stimdata.muscle(i).target(j).rep(k).meanWristPosition_global];
            stimdata.muscle(i).target(j).allMeanWristQuaternion     = [stimdata.muscle(i).target(j).allMeanWristQuaternion; stimdata.muscle(i).target(j).rep(k).meanWristQuaternion];
            
            % find mean force in thorax frame
            trot=qGetR(stimdata.muscle(i).target(j).rep(k).meanThoraxQuaternion);
            f= stimdata.muscle(i).target(j).rep(k).meanNewHmForce;
            force_thor=inv(trot)*f';
            force_thor=force_thor';
            stimdata.muscle(i).target(j).allMeanNewHmForce_Thor     = [stimdata.muscle(i).target(j).allMeanNewHmForce_Thor; force_thor];
        end
        stimdata.muscle(i).target(j).meanThoraxQuaternion       = nanmean(stimdata.muscle(i).target(j).allThoraxQuaternion);
        stimdata.muscle(i).target(j).meanHmMeasPosition         = nanmean(stimdata.muscle(i).target(j).allHmMeasPosition);
        stimdata.muscle(i).target(j).meanWristPosition          = nanmean(stimdata.muscle(i).target(j).allWristPosition);
        stimdata.muscle(i).target(j).meanWristPosition_global   = nanmean(stimdata.muscle(i).target(j).allWristPosition_global);
        stimdata.muscle(i).target(j).meanWristQuaternion          = nanmean(stimdata.muscle(i).target(j).allWristQuaternion);
        
        
    end
    
    
end


save('parseddata.mat','stimdata')