%function [angles,velocities,accelerations,torqueFilt,newHumerusPosition,newThoraxRotation,humerusStartIndex,humerusEndIndex,jointAngles,masMarkers] = smoothmotion(time,forceFilt,thoraxPosition,thoraxRotation,thoraxFlag,humerusPosition,humerusRotation,humerusFlag,forearmPosition,forearmRotation,forearmFlag,hmPosition,hmRotation,hmFlag,newHmMeasPosition,masMarkers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHMOTION_Wrist.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: This function calculates the wrist position relative to the
%   thorax frame. Moving average filter is used to remove noise.
%       Built off SMOOTHMOTION.M
%       Only 2 rigid bodies. Forearm and Thorax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
% Ouputs:
%   filteredWristPosition
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author of SMOOTHMOTION: Eric Schearer
% Author of SMOOTHMOTION_Wrist: Derek Wolf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 03/09/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% constant parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
plots = 1;


% start at the first visible forearm flag and end at the last one
thoraxStartIndex    = find(thoraxFlag==1,1,'first');
forearmStartIndex   = find(forearmFlag==1,1,'first');
bigStartIndex       = find(thoraxFlag==1 & forearmFlag==1,1,'first');
bigEndIndex         = find(thoraxFlag==1 & forearmFlag==1,1,'last');

forearmEndIndex = find(forearmFlag==1,1,'last');

% if thorax is never visible. Set first element of thorax equal to thoraxQ
if isempty(thoraxStartIndex)
    thoraxStartIndex=5;
    thoraxQuaternion(thoraxStartIndex,:)=thoraxQ;
    thoraxPosition(thoraxStartIndex,:)=thoraxP;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check rigid body flags and replace data with nan if flag = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(time)
    if thoraxFlag(i) == 0
        thoraxPosition(i,:) = thoraxPosition(thoraxStartIndex,:);
        thoraxQuaternion(i,:) = thoraxQuaternion(thoraxStartIndex,:);
        thoraxAxisAngle(i,:) = quat2axisAngle(thoraxQuaternion(thoraxStartIndex,:));
        jointAngles(i,1:3) = NaN*ones(1,3);
    else
        thoraxAxisAngle(i,:) = quat2axisAngle(thoraxQuaternion(i,:));
    end
    if forearmFlag(i) == 0
        forearmPosition(i,:) = NaN*ones(1,3);
        forearmQuaternion(i,:) = NaN*ones(1,4);
        forearmAxisAngle(i,:) = NaN*ones(1,3);
        jointAngles(i,4:6) = NaN*ones(1,3);
    else
        forearmAxisAngle(i,:) = quat2axisAngle(forearmQuaternion(i,:));
    end
end

% Get pose of forearm in thorax frame
for i=1:length(time)
    % Position of forearm in thorax frame
    pos=[forearmPosition(i,:)';1];
    Tr=qGetR(thoraxQuaternion(i,:));
    Tp=thoraxPosition(i,:)';
    invT=[inv(Tr) -inv(Tr)*Tp];
    invT=[invT;0 0 0 1];
    newpos=invT*pos;
    forearmPos_thor(i,:)=newpos(1:3);
    
    % Keeps forearmPosition in global frame
    forearmPos_global(i,:)=forearmPosition(i,:);
    
    % Orientation of forearm in thorax frame
    frot=qGetR(forearmQuaternion(i,:));
    trot=qGetR(thoraxQuaternion(i,:));
    forearmQuaternion_thor(i,:)=forearmQuaternion(i,:);
end




filteredWristPosition = [ones(forearmStartIndex-1,1)*forearmPos_thor(1,:);
    forearmPos_thor(forearmStartIndex:forearmEndIndex,:);
    ones(length(time)-forearmEndIndex,1)*forearmPos_thor(end,:)];

filteredWristPosition_global = [ones(forearmStartIndex-1,1)*forearmPos_thor(1,:);
    forearmPos_thor(forearmStartIndex:forearmEndIndex,:);
    ones(length(time)-forearmEndIndex,1)*forearmPos_global(end,:)];

if plots == 1
    % Plot wrist position
    figure(1);
    plot(filteredWristPosition);   
end


