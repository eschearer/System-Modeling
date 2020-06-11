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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%
% define time parameters
%%%%%%%%%%%%%%%%%%%%%%%%
stimTime = 1.25;
restTime = 0.75;
riseTime = 0.25;
global dt;
dt = 1/52;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the data structure for storing things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muscles = 10;
targets = 27;
reps = 5;
for i = 1:muscles
    for j = 1:targets
        for k = 1:reps
            stimdata.muscle(i).target(j).rep(k).jointAngles         = zeros(1/dt:6);
            stimdata.muscle(i).target(j).rep(k).filteredJointAngles = zeros(1/dt:6);
            stimdata.muscle(i).target(j).rep(k).handMarker          = zeros(1/dt:4);
            stimdata.muscle(i).target(j).rep(k).hmForce             = zeros(1/dt:3);
            stimdata.muscle(i).target(j).rep(k).newHmForce          = zeros(1/dt:3);
            stimdata.muscle(i).target(j).rep(k).hmMeasPosition      = zeros(1/dt:3);
            stimdata.muscle(i).target(j).rep(k).stimLevel           = zeros(1/dt:1);
            stimdata.muscle(i).target(j).rep(k).stimPW              = zeros(1/dt:12);
            stimdata.muscle(i).target(j).rep(k).thoraxQuaternion    = zeros(1/dt:4);
            stimdata.muscle(i).target(j).rep(k).jointTorques        = zeros(1/dt:5);
            stimdata.muscle(i).target(j).rep(k).wristPosition       = zeros(1/dt:3);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% load the target order
%%%%%%%%%%%%%%%%%%%%%%%
load('../../executionCode/targetOrder1.txt')
load('../../executionCode/targetOrder2.txt')
load('../../executionCode/targetOrder3.txt')
targetOrder = [targetOrder1 targetOrder2 targetOrder3];
repOrder = [ones(targets,1) 2*ones(targets,1) 3*ones(targets,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through all of the data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trials = 81;
trialIndex = [1:81];
for l = 1:length(trialIndex)
    %close all;
    trial = trialIndex(l);
    load(['dataTrial',num2str(trial)])
    trial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % use RTS smoother to find joint angles 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    smoothmotion
    pause 
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
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).jointAngles        = jointAngles(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).filteredJointAngles= filteredJointAngles(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).handMarker         = handMarker(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).hmForce            = hmForce(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).newHmForce         = newHmForce(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).hmMeasPosition     = hmMeasPosition(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).stimLevel          = stimLevel(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).stimPW             = stimPW(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).thoraxQuaternion   = thoraxQuaternion(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:);
        stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).wristPosition      = filteredWristPosition(startIndex(i)+riseTime/dt:startIndex(i)+stimTime/dt,:); 
%         muscle(startIndex(i))
%         figure(1)
%         hold on
%         plot(stimdata.muscle(muscle(startIndex(i))).target(currentTarget).rep(currentRep).hmForce(:,1))
        
        
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
        stimdata.muscle(i).target(j).allJointAngles         = [];
        stimdata.muscle(i).target(j).allFilteredJointAngles = [];
        stimdata.muscle(i).target(j).allHandMarker          = [];
        stimdata.muscle(i).target(j).allHmForce             = [];
        stimdata.muscle(i).target(j).allNewHmForce          = [];
        stimdata.muscle(i).target(j).allHmMeasPosition      = [];
        stimdata.muscle(i).target(j).allStimLevel           = [];
        stimdata.muscle(i).target(j).allStimPW              = [];
        stimdata.muscle(i).target(j).allThoraxQuaternion    = [];
        stimdata.muscle(i).target(j).allJointTorques        = [];
        stimdata.muscle(i).target(j).allWristPosition       = [];
        for k = 1:reps
            stimdata.muscle(i).target(j).allJointAngles         = [stimdata.muscle(i).target(j).allJointAngles; stimdata.muscle(i).target(j).rep(k).jointAngles];
            stimdata.muscle(i).target(j).allFilteredJointAngles = [stimdata.muscle(i).target(j).allFilteredJointAngles; stimdata.muscle(i).target(j).rep(k).filteredJointAngles];
            stimdata.muscle(i).target(j).allHandMarker          = [stimdata.muscle(i).target(j).allHandMarker; stimdata.muscle(i).target(j).rep(k).handMarker];
            stimdata.muscle(i).target(j).allHmForce             = [stimdata.muscle(i).target(j).allHmForce; stimdata.muscle(i).target(j).rep(k).hmForce];
            stimdata.muscle(i).target(j).allNewHmForce          = [stimdata.muscle(i).target(j).allNewHmForce; stimdata.muscle(i).target(j).rep(k).newHmForce];
            stimdata.muscle(i).target(j).allHmMeasPosition      = [stimdata.muscle(i).target(j).allHmMeasPosition; stimdata.muscle(i).target(j).rep(k).hmMeasPosition];
            stimdata.muscle(i).target(j).allStimLevel           = [stimdata.muscle(i).target(j).allStimLevel; stimdata.muscle(i).target(j).rep(k).stimLevel];
            stimdata.muscle(i).target(j).allStimPW              = [stimdata.muscle(i).target(j).allStimPW; stimdata.muscle(i).target(j).rep(k).stimPW];
            stimdata.muscle(i).target(j).allThoraxQuaternion    = [stimdata.muscle(i).target(j).allThoraxQuaternion; stimdata.muscle(i).target(j).rep(k).thoraxQuaternion];
            stimdata.muscle(i).target(j).allWristPosition       = [stimdata.muscle(i).target(j).allWristPosition; stimdata.muscle(i).target(j).rep(k).wristPosition];
        end
        stimdata.muscle(i).target(j).meanJointAngles            = nanmean(stimdata.muscle(i).target(j).allJointAngles);
        stimdata.muscle(i).target(j).meanFilteredJointAngles    = nanmean(stimdata.muscle(i).target(j).allFilteredJointAngles);
        stimdata.muscle(i).target(j).meanThoraxQuaternion       = nanmean(stimdata.muscle(i).target(j).allThoraxQuaternion);
        stimdata.muscle(i).target(j).meanHmMeasPosition         = nanmean(stimdata.muscle(i).target(j).allHmMeasPosition);
        stimdata.muscle(i).target(j).meanWristPosition          = nanmean(stimdata.muscle(i).target(j).allWristPosition);
       
        
        % compute torques if you have data to do so
        if isnan(stimdata.muscle(i).target(j).meanFilteredJointAngles) ~= 1
            thoraxRotMat        = qGetR(stimdata.muscle(i).target(j).meanThoraxQuaternion);
            ELEVATIONPLANE      = stimdata.muscle(i).target(j).meanFilteredJointAngles(1);
            SHOULDERELEVATION   = stimdata.muscle(i).target(j).meanFilteredJointAngles(2);
            SHOULDERROTATION    = stimdata.muscle(i).target(j).meanFilteredJointAngles(3);
            ELBOWFLEXION        = stimdata.muscle(i).target(j).meanFilteredJointAngles(4);
            ELBOWCARRY          = stimdata.muscle(i).target(j).meanFilteredJointAngles(5);
            ELBOWPRONATION      = stimdata.muscle(i).target(j).meanFilteredJointAngles(6);

            ballPos = stimdata.muscle(i).target(j).meanHmMeasPosition; % position of the ball and socket joint
            wristPos = stimdata.muscle(i).target(j).meanWristPosition;
            jac = compute_Jall(ELBOWCARRY,ELBOWFLEXION,ELEVATIONPLANE,FOREARMLENGTH,HUMERUSLENGTH,SHOULDERROTATION,SHOULDERELEVATION,thoraxRotMat(1,1),thoraxRotMat(1,2),thoraxRotMat(1,3),thoraxRotMat(2,1),thoraxRotMat(2,2),thoraxRotMat(2,3),thoraxRotMat(3,1),thoraxRotMat(3,2),thoraxRotMat(3,3));

            % loop through all the torques for the muscle and target
            for count = 1:length(stimdata.muscle(i).target(j).allWristPosition);
                extraTorque = cross(ballPos'-wristPos',stimdata.muscle(i).target(j).allHmForce(count,:)); % N-mm resultant torque at the wrist given a force is applied at the ball and socket joint
                stimdata.muscle(i).target(j).allJointTorques(count,:) = (jac'*[stimdata.muscle(i).target(j).allHmForce(count,:)';extraTorque'])'/1000; % N-m     
            end
        end
    end
    
    
end


save('parseddata.mat','stimdata')