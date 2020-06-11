%function [angles,velocities,accelerations,torqueFilt,newHumerusPosition,newThoraxRotation,humerusStartIndex,humerusEndIndex,jointAngles,masMarkers] = smoothmotion(time,forceFilt,thoraxPosition,thoraxRotation,thoraxFlag,humerusPosition,humerusRotation,humerusFlag,forearmPosition,forearmRotation,forearmFlag,hmPosition,hmRotation,hmFlag,newHmMeasPosition,masMarkers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHMOTION.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Descritpion: This function uses an RTS smoother to smooth motion capture
% data.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: 
%   
% Ouputs:
%   filteredJointAngles
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 8 November 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start at the first visable humerus flag and end at the last one 
humerusStartIndex   = find(humerusFlag==1,1,'first');
thoraxStartIndex    = find(thoraxFlag==1,1,'first');
forearmStartIndex   = find(forearmFlag==1,1,'first');
bigStartIndex       = find(humerusFlag==1 & forearmFlag==1,1,'first');
bigEndIndex         = find(humerusFlag==1 & forearmFlag==1,1,'last');
% if humerusStartIndex > 0.5/dt
%     humerusStartIndex = humerusStartIndex - 0.5/dt;
% end
humerusEndIndex = find(humerusFlag==1,1,'last');
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
        % recompute the joint angles after with the good thorax orientation
        %if humerusFlag(i)
        %    thoraxRotTemp = qGetR(thoraxQuaternion(i,:));
        %    humerusRotTemp = qGetR(humerusQuaternion(i,:));
            %[jointAngles(i,1), jointAngles(i,2), jointAngles(i,3)] = rotyzy(thoraxRotTemp'*humerusRotTemp);
        %end
    end
    if humerusFlag(i) == 0
        humerusPosition(i,:) = NaN*ones(1,3);
        humerusQuaternion(i,:) = NaN*ones(1,4);
        humerusAxisAngle(i,:) = NaN*ones(1,3);
        jointAngles(i,1:6) = NaN*ones(1,6);
    else
        humerusAxisAngle(i,:) = quat2axisAngle(humerusQuaternion(i,:));
    end
    if forearmFlag(i) == 0
        forearmPosition(i,:) = NaN*ones(1,3);
        forearmQuaternion(i,:) = NaN*ones(1,4);
        forearmAxisAngle(i,:) = NaN*ones(1,3);
        jointAngles(i,4:6) = NaN*ones(1,3);
    else
        forearmAxisAngle(i,:) = quat2axisAngle(forearmQuaternion(i,:));
    end
    if handMarker(i,4) == 0
        handMarker(i,1:3) = NaN*ones(1,3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% constant parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
FOREARMLENGTH = 269;
HUMERUSLENGTH = 306;
ELEMdist = 59;
RSUSdist = 55;
plots = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kinematic tree
% parent jointType sensorType jointAxis jointCategory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
segment = [...
   0 0 0 0 0;  %1   root
   1 2 0 2 1;  %2   thorax position
   2 4 0 2 1;  %3   thorax orientation
   3 2 0 0 1;  %4   shoulder position
   4 3 0 2 1;  %5   shoulder elevation plane
   5 3 0 3 1;  %6   shoulder elevation
   6 3 0 2 1;  %7   shoulder rotation
   7 1 0 2 2;  %8   humerus length
   8 3 0 1 1;  %9   elbow flexion
   9 3 0 3 2;  %10  elbow carry
  10 3 0 2 1;  %11  elbow pronation
  11 1 0 2 2;  %12  forearm length
  %12 2 0 0 2;  %13  hand marker        
   2 2 2 0 2;  %14  pos thorax
   3 4 4 0 2;  %15  quat thorax
   4 2 2 0 2;  %16  pos humerus
   7 4 4 0 2;  %17  quat humerus
   8 2 2 0 2;  %18  pos forearm
  11 4 4 0 2;  %19  quat forearm
  %13 2 2 0 2;  %20  Hand marker position
];

% noise and dynamics model (compressed)
ratio = 10;
global dt
dt = 1/52;
datR = [0, 20];
datS = [3, 20];
datV = 0.2;
mode = 1;       % estimate both variable and constant state parameters

[map, info, S, R, V] = prepareJerk(segment, mode, dt,datS,datR,datV,ratio);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the output vector for the optimal state estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = [thoraxPosition(humerusStartIndex:humerusEndIndex,:)'; 
     thoraxQuaternion(humerusStartIndex:humerusEndIndex,:)'; 
     humerusPosition(humerusStartIndex:humerusEndIndex,:)'; 
     humerusQuaternion(humerusStartIndex:humerusEndIndex,:)';
     forearmPosition(humerusStartIndex:humerusEndIndex,:)'; 
     forearmQuaternion(humerusStartIndex:humerusEndIndex,:)'...
     %handMarker(humerusStartIndex:humerusEndIndex,1:3)'...
     ];
keepTime = time(humerusStartIndex:humerusEndIndex);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build the sensor covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for opto and HM sensors
thoraxPosVar = 10;%100;
thoraxRotVar = 0.1;%0.5;
humerusPosVar = 100;%100;
humerusRotVar = 1;%0.5;
forearmPosVar = 100;%100;
forearmRotVar = 1;%0.5;
markerVariance = 10;%100;
hmVariance = 1e-3;%;%10;

V = diag([thoraxPosVar*ones(1,3)...
          thoraxRotVar*ones(1,3)...
          humerusPosVar*ones(1,3)... 
          humerusRotVar*ones(1,3)...
          forearmPosVar*ones(1,3)... 
          forearmRotVar*ones(1,3)...
          %markerVariance*ones(1,3)
          ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build the process covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smallProcVar = 1e-2;
bigProcVar = 1e-1;
thoraxPosProcVar        = smallProcVar;
thoraxRotProcVar        = smallProcVar;
humerusPosProcVar       = smallProcVar;
humerusRotProcVar       = smallProcVar;
humerusLengthProcVar    = smallProcVar;
forearmRotProcVar       = smallProcVar;
forearmLengthProcVar    = smallProcVar;
hmPosProcVar            = smallProcVar;
sensPosProcVar          = smallProcVar;
thoraxVelProcVar        = smallProcVar;
thoraxRotVelProcVar     = smallProcVar;
humerusVelProcVar       = smallProcVar;
humerusRotVelProcVar    = smallProcVar;
humerusLengthVelProcVar = smallProcVar;
forearmRotVelProcVar    = smallProcVar;
hmVelProcVar            = smallProcVar;
thoraxAccProcVar        = smallProcVar;
thoraxRotAccProcVar     = smallProcVar;
humerusAccProcVar       = smallProcVar;
humerusRotAccProcVar    = smallProcVar;
humerusLengthAccProcVar = smallProcVar;
forearmRotAccProcVar    = smallProcVar;
hmAccProcVar            = smallProcVar;
thoraxJerkProcVar       = bigProcVar;
thoraxRotJerkProcVar    = bigProcVar;
humerusJerkProcVar      = bigProcVar;
humerusRotJerkProcVar   = bigProcVar;
humerusLengthJerkProcVar = smallProcVar;
forearmRotJerkProcVar   = bigProcVar;
hmJerkProcVar           = bigProcVar;
hmSensPosProcVar        = smallProcVar;

%R = diag([thoraxPosProcVar*ones(1,3) thoraxRotProcVar*ones(1,3) humerusPosProcVar*ones(1,3) humerusRotProcVar*ones(1,3) humerusLengthProcVar forearmRotProcVar*ones(1,3) forearmLengthProcVar hmPosProcVar*ones(1,3) sensPosProcVar*ones(1,18) hmSensPosProcVar*ones(1,3) humerusVelProcVar*ones(1,3) humerusRotVelProcVar*ones(1,3) forearmRotVelProcVar*ones(1,2) humerusAccProcVar*ones(1,3) humerusRotAccProcVar*ones(1,3) forearmRotAccProcVar*ones(1,2) humerusJerkProcVar*ones(1,3) humerusRotJerkProcVar*ones(1,3) forearmRotJerkProcVar*ones(1,2)]);
R = diag([thoraxPosProcVar*ones(1,3)...
          thoraxRotProcVar*ones(1,3)...
          humerusPosProcVar*ones(1,3)... 
          humerusRotProcVar*ones(1,3)... 
          humerusLengthProcVar... 
          forearmRotProcVar*ones(1,3)... 
          forearmLengthProcVar...%hmPosProcVar*ones(1,3)...
          sensPosProcVar*ones(1,18)...
          thoraxVelProcVar*ones(1,3)...
          thoraxRotVelProcVar*ones(1,3)...
          humerusVelProcVar*ones(1,3)...
          humerusRotVelProcVar*ones(1,3)... 
          forearmRotVelProcVar*ones(1,2)...
          thoraxAccProcVar*ones(1,3)...
          thoraxRotAccProcVar*ones(1,3)...
          humerusAccProcVar*ones(1,3)... 
          humerusRotAccProcVar*ones(1,3)...
          forearmRotAccProcVar*ones(1,2)...
          thoraxJerkProcVar*ones(1,3)...
          thoraxRotJerkProcVar*ones(1,3)...
          humerusJerkProcVar*ones(1,3)...
          humerusRotJerkProcVar*ones(1,3)...
          forearmRotJerkProcVar*ones(1,2)...
          ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build the initial state covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conAngleCov = 0.001;
varAngleCov = 10000;
conPosCov = 0.1;
varPosCov = 1000000;
posCov = 400;
angleCov = 1000;% 0.1;
velCov = 0.001;%1;
angVelCov = 0.0001; 
accCov = 0.0001;%0.1;
angAccCov = 0.00001; 
jerkCov = 0.00001;%0.1;
angJerkCov = 0.00001; 

thoraxPosCov = conPosCov*ones(1,3);
thoraxRotCov = conAngleCov*ones(1,3);
humerusPosCov = varPosCov*ones(1,3);
humerusRotCov = varAngleCov*ones(1,3);
humerusLengthCov = conPosCov;
forearmRotCov = varAngleCov*ones(1,3);
forearmLengthCov = conPosCov;
handMarkerPosCov = conPosCov*ones(1,3);
sensorPositionCov = conPosCov*ones(1,18);
thoraxVelCov = velCov*ones(1,3);
thoraxRotVelCov = angVelCov*ones(1,3);
humerusVelCov = velCov*ones(1,3);
humerusRotVelCov = angVelCov*ones(1,3);
humerusLengthVelCov = conPosCov;
forearmRotVelCov = angVelCov*ones(1,2);
handMarkerVelCov = conPosCov*ones(1,3);
thoraxAccCov = accCov*ones(1,3);
thoraxRotAccCov = angAccCov*ones(1,3);
humerusAccCov = accCov*ones(1,3);
humerusRotAccCov = angAccCov*ones(1,3);
humerusLengthAccCov = conPosCov;
forearmRotAccCov = angAccCov*ones(1,2);
handMarkerAccCov = conPosCov*ones(1,3);
thoraxJerkCov = jerkCov*ones(1,3);
thoraxRotJerkCov = angJerkCov*ones(1,3);
humerusJerkCov = jerkCov*ones(1,3);
humerusRotJerkCov = angJerkCov*ones(1,3);
humerusLengthJerkCov = conPosCov;
forearmRotJerkCov = angJerkCov*ones(1,2);
handMarkerJerkCov = conPosCov*ones(1,3);

S = diag([thoraxPosCov...
          thoraxRotCov...
          humerusPosCov...
          humerusRotCov...
          humerusLengthCov... 
          forearmRotCov... 
          forearmLengthCov... %handMarkerPosCov...
          sensorPositionCov... 
          thoraxVelCov... 
          thoraxRotVelCov...
          humerusVelCov... 
          humerusRotVelCov...
          forearmRotVelCov...
          thoraxAccCov... 
          thoraxRotAccCov...
          humerusAccCov... 
          humerusRotAccCov... 
          forearmRotAccCov...
          thoraxJerkCov... 
          thoraxRotJerkCov...
          humerusJerkCov... 
          humerusRotJerkCov...
          forearmRotJerkCov...
          ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thoraxPos0              = thoraxPosition(thoraxStartIndex,:);
thoraxRot0              = quat2axisAngle(thoraxQuaternion(thoraxStartIndex,:)); 
thoraxRotMat            = qGetR(thoraxQuaternion(thoraxStartIndex,:));
% find the first visible humerus flag and use this data as the
% initial guess
humerusPos0             = thoraxRotMat'*(humerusPosition(humerusStartIndex,:)-thoraxPosition(thoraxStartIndex,:))';
% find the first non-NaN value for the shoulder joint angles
shoulderIndex           = find(~isnan(jointAngles(:,1)),1,'first');%find(humerusFlag==1.*thoraxFlag==1,1,'first');
elevationPlane0         =  jointAngles(shoulderIndex,1)*pi/180;
shoulderElevation0      =  jointAngles(shoulderIndex,2)*pi/180;
shoulderRotation0       =  jointAngles(shoulderIndex,3)*pi/180;
humerusLength0          = -HUMERUSLENGTH;
% find the first non-NaN value for the elbow joint angles
elbowIndex              = find(~isnan(jointAngles(:,4)),1,'first');%find(humerusFlag==1.*forearmFlag==1,1,'first');
elbowFlexion0           = jointAngles(elbowIndex,4)*pi/180;
elbowCarry0             = nanmean(jointAngles(:,5))*pi/180;
elbowPronation0         = jointAngles(elbowIndex,6)*pi/180;
forearmLength0          = -FOREARMLENGTH;
handMarker0             = [-10 0 -10];
hmPos0                  = [0 -80 -60];%[0 50 -70];
hmVel0                  = [0 0 0];
hmAcc0                  = [0 0 0];
% assume that thorax, humerus, and forearm optotrak rigid bodies were
% correctly positioned and aligned
thoraxPosSens0          = [0 0 0];
thoraxRotSens0          = [0 0 0];
humerusPosSens0         = [0 0 0];
humerusRotSens0         = [0 0 0];
forearmPosSens0         = [0 0 0];
forearmRotSens0         = [0 0 0];
thoraxVelocity0         = [0 0 0];
thoraxRotVelocity0      = [0 0 0];
humerusVelocity0        = [0 0 0];
elevationPlaneVelocity0 = 0;
shoulderElevationVelocity0  = 0;
shoulderRotationVelocity0   = 0;
humerusLengthVelocity0      = 0;
elbowFlexionVelocity0       = 0;
elbowPronationVelocity0     = 0;
handMarkerVel0             = [0 0 0];
thoraxAcceleration0         = [0 0 0];
thoraxRotAccceleration0     = [0 0 0];
humerusAcceleration0        = [0 0 0];
elevationPlaneAcceleration0 = 0;
shoulderElevationAcceleration0  = 0;
shoulderRotationAcceleration0   = 0;
humerusLengthAcceleration0      = 0;
elbowFlexionAcceleration0       = 0;
elbowPronationAcceleration0     = 0;
handMarkerAcc0             = [0 0 0];
thoraxJerk0                     = [0 0 0];
thoraxRotJerk0                  = [0 0 0];
humerusJerk0                    = [0 0 0];
elevationPlaneJerk0             = 0;
shoulderElevationJerk0          = 0;
shoulderRotationJerk0           = 0;
humerusLengthJerk0              = 0;
elbowFlexionJerk0               = 0;
elbowPronationJerk0             = 0;
handMarkerJerk0             = [0 0 0];
GHposSens0                      = [0 0 0];
ELposSens0                      = [ELEMdist/2 0 0];
EMposSens0                      = [-ELEMdist/2 0 0];
RSposSens0                      = [RSUSdist/2 0 0];
USposSens0                      = [-RSUSdist/2 0 0];
hmPosSens0                      = [0 0 0];%

% stack the initial states into a vector
x0 = [...
    thoraxPos0';thoraxRot0';humerusPos0;
    elevationPlane0;shoulderElevation0;shoulderRotation0;humerusLength0;
    elbowFlexion0;elbowCarry0;elbowPronation0;forearmLength0;
    %handMarker0';
    thoraxPosSens0';thoraxRotSens0';
    humerusPosSens0';humerusRotSens0';
    forearmPosSens0';forearmRotSens0';
    %hmPosSens0';
    thoraxVelocity0';
    thoraxRotVelocity0';
    humerusVelocity0';
    elevationPlaneVelocity0;shoulderElevationVelocity0;shoulderRotationVelocity0;
    %humerusLengthVelocity0;
    elbowFlexionVelocity0;elbowPronationVelocity0;
    thoraxAcceleration0';
    thoraxRotAccceleration0';
    humerusAcceleration0';
    elevationPlaneAcceleration0;shoulderElevationAcceleration0;shoulderRotationAcceleration0;
    %humerusLengthAcceleration0;
    elbowFlexionAcceleration0;elbowPronationAcceleration0;
    thoraxJerk0';
    thoraxRotJerk0';
    humerusJerk0';
    elevationPlaneJerk0;shoulderElevationJerk0;shoulderRotationJerk0;
    %humerusLengthJerk0;
    elbowFlexionJerk0;elbowPronationJerk0...
    ];
  
nS = size(segment,1);
A = diag(ones(map.nX,1));
for i = 1:nS
    if map.aP(i) > 0 
        switch segment(i,2)
            case 1 % sliding joint
                A(map.aP(i)+map.aW(i),map.aD(i)) = dt;
                A(map.aD(i),map.aDD(i)) = dt;
                A(map.aDD(i),map.aDDD(i)) = dt;
                A(map.aDDD(i),map.aDDD(i)) = 0;
            case 2 % translation joint
                A(map.aP(i)+map.aW(i),map.aD(i)) = dt;
                A(map.aP(i)+map.aW(i)+1,map.aD(i)+1) = dt;
                A(map.aP(i)+map.aW(i)+2,map.aD(i)+2) = dt;
                A(map.aD(i),map.aDD(i)) = dt;
                A(map.aD(i)+1,map.aDD(i)+1) = dt;
                A(map.aD(i)+2,map.aDD(i)+2) = dt;
                A(map.aDD(i),map.aDDD(i)) = dt;
                A(map.aDD(i)+1,map.aDDD(i)+1) = dt;
                A(map.aDD(i)+2,map.aDDD(i)+2) = dt;
                A(map.aDDD(i),map.aDDD(i)) = 0;
                A(map.aDDD(i)+1,map.aDDD(i)+1) = 0;
                A(map.aDDD(i)+2,map.aDDD(i)+2) = 0;
            case 3 % hinge joint
                A(map.aP(i)+map.aW(i),map.aD(i)) = dt;
                A(map.aD(i),map.aDD(i)) = dt;
                A(map.aDD(i),map.aDDD(i)) = dt;
                A(map.aDDD(i),map.aDDD(i)) = 0;
            case 4 % rotation joint
                A(map.aP(i)+map.aW(i),map.aD(i)) = dt;
                A(map.aP(i)+map.aW(i)+1,map.aD(i)+1) = dt;
                A(map.aP(i)+map.aW(i)+2,map.aD(i)+2) = dt;
                A(map.aD(i),map.aDD(i)) = dt;
                A(map.aD(i)+1,map.aDD(i)+1) = dt;
                A(map.aD(i)+2,map.aDD(i)+2) = dt;
                A(map.aDD(i),map.aDDD(i)) = dt;
                A(map.aDD(i)+1,map.aDDD(i)+1) = dt;
                A(map.aDD(i)+2,map.aDDD(i)+2) = dt;
                A(map.aDDD(i),map.aDDD(i)) = 0;
                A(map.aDDD(i)+1,map.aDDD(i)+1) = 0;
                A(map.aDDD(i)+2,map.aDDD(i)+2) = 0;
        end
    end
end
%%
% do the forward extended Kalman filter
[X, Xminus, Err, S, Pminus, Pplus, Res, Pred, Cost] = estimate_ekf(Y, x0, S, segment, map, info, R, V, A);
% do the RTS smoother
[Xsmooth, Psmooth] = rtssmoother(Xminus,X,Pminus,Pplus,A);

elbow = zeros(length(Xsmooth),3);
wrist = zeros(length(Xsmooth),3);
hm = zeros(length(Xsmooth),3);
for i = 1:length(Xsmooth)
    [frmPos, frmQuat, frmRot, J] = kinematics(Y(:,i), Xsmooth(:,i), segment, map); 
    hm(i,:) = frmPos(:,13)';
    wrist(i,:) = frmPos(:,12)';
end

if plots == 2
    
%     figure(1)
%     hold on
%     plot(time,thoraxPosition)
%     plot(keepTime,Xsmooth(1:3,:),'LineWidth',4)
    
    figure(1)
    hold on
    plot(keepTime,jointAngles(humerusStartIndex:humerusEndIndex,1),'b')
    plot(keepTime,jointAngles(humerusStartIndex:humerusEndIndex,2),'g')
    plot(keepTime,jointAngles(humerusStartIndex:humerusEndIndex,3),'r')
    plot(keepTime,Xsmooth(10,:)*180/pi','b','LineWidth',4)
    plot(keepTime,Xsmooth(11,:)*180/pi','g','LineWidth',4)
    plot(keepTime,Xsmooth(12,:)*180/pi','r','LineWidth',4)
    title('shoulder joint angles')

    figure(2)
    hold on
    plot(keepTime,Xsmooth(map.aD(5),:)','b','LineWidth',4)
    plot(keepTime,Xsmooth(map.aD(6),:)','g','LineWidth',4)
    plot(keepTime,Xsmooth(map.aD(7),:)','r','LineWidth',4)
    plot(keepTime,[0 diff(Xsmooth(10,:))/dt],'b','LineWidth',2)
    plot(keepTime,[0 diff(Xsmooth(11,:))/dt],'g','LineWidth',2)
    plot(keepTime,[0 diff(Xsmooth(12,:))/dt],'r','LineWidth',2)
    title('shoulder joint velocities')

%     figure(3)
%     hold on
%     plot(keepTime,Xsmooth(map.aDD(5),:)','b','LineWidth',4)
%     plot(keepTime,Xsmooth(map.aDD(6),:)','g','LineWidth',4)
%     plot(keepTime,Xsmooth(map.aDD(7),:)','r','LineWidth',4)
%     plot(keepTime,[0 diff(Xsmooth(map.aD(5),:))/dt],'b','LineWidth',2)
%     plot(keepTime,[0 diff(Xsmooth(map.aD(6),:))/dt],'g','LineWidth',2)
%     plot(keepTime,[0 diff(Xsmooth(map.aD(7),:))/dt],'r','LineWidth',2)
%     title('shoulder joint accelerations')
% 
    figure(4)
    hold on
    plot(keepTime,jointAngles(humerusStartIndex:humerusEndIndex,4)','b')
    plot(keepTime,jointAngles(humerusStartIndex:humerusEndIndex,5)','g')
    plot(keepTime,jointAngles(humerusStartIndex:humerusEndIndex,6)','r')
    plot(keepTime,Xsmooth(14,:)*180/pi','b','LineWidth',4)
    plot(keepTime,Xsmooth(15,:)*180/pi','g','LineWidth',4)
    plot(keepTime,Xsmooth(16,:)*180/pi','r','LineWidth',4)
    title('elbow joint angles')
% 
    figure(5)
    hold on
    plot(keepTime,Xsmooth(map.aD(9),:)','b','LineWidth',4)
    plot(keepTime,Xsmooth(map.aD(11),:)','r','LineWidth',4)
    plot(keepTime,[0 diff(Xsmooth(14,:))/dt],'b','LineWidth',2)
    plot(keepTime,[0 diff(Xsmooth(16,:))/dt],'r','LineWidth',2)    
    title('elbow joint velocities')
% 
%     figure(6)
%     hold on
%     plot(keepTime,Xsmooth(map.aDD(9),:)','b','LineWidth',4)
%     plot(keepTime,Xsmooth(map.aDD(11),:)','r','LineWidth',4)
%     plot(keepTime,[0 diff(Xsmooth(map.aD(9),:))/dt],'b','LineWidth',2)
%     plot(keepTime,[0 diff(Xsmooth(map.aD(11),:))/dt],'r','LineWidth',2)
%     title('elbow joint accelerations')
end
% if plots == 2
%     figure(14)
%     hold on
%     plot(keepTime,hm(:,1),'r','LineWidth',3)
%     plot(keepTime,hm(:,2),'g','LineWidth',3)
%     plot(keepTime,hm(:,3),'b','LineWidth',3)
%     plot(keepTime,handMarker(humerusStartIndex:humerusEndIndex,1),'r','LineWidth',1)
%     plot(keepTime,handMarker(humerusStartIndex:humerusEndIndex,2),'g','LineWidth',1)
%     plot(keepTime,handMarker(humerusStartIndex:humerusEndIndex,3),'b','LineWidth',1)
%     title('Hand Position')
%     
%     
% end

% rep1.torqueFilt = zeros(size(Xsmooth,2),4);
% rep2.torqueFilt = zeros(size(Xsmooth,2),4);
% %rep3.torqueFilt = zeros(size(Xsmooth,2),4);
% rep1.extraTorque = zeros(size(Xsmooth,2),3);
% 
% saveJall = zeros(6,5,size(Xsmooth,2));
% for i = 1:size(Xsmooth,2)
%     ballPos = hmMeasPosition(humerusStartIndex+1-1,:); % position of the ball and socket joint
%     wristPos = frmPos(:,12);
%     thoraxRot = thoraxRotMat; % qGetR(thoraxQuaternion(i,:));
%     Jall = compute_Jall(Xsmooth(15,i),Xsmooth(14,i),Xsmooth(10,i),FOREARMLENGTH,HUMERUSLENGTH,Xsmooth(12,i),Xsmooth(11,i),thoraxRot(1,1),thoraxRot(1,2),thoraxRot(1,3),thoraxRot(2,1),thoraxRot(2,2),thoraxRot(2,3),thoraxRot(3,1),thoraxRot(3,2),thoraxRot(3,3)); % mm
%     saveJall(:,:,i) = Jall;
%     rep1.extraTorque = cross(ballPos'-wristPos,rep1.forceFilt(humerusStartIndex+i-1,:)); % N-mm resultant torque at the wrist given a force is applied at the ball and socket joint
%     rep2.extraTorque = cross(ballPos'-wristPos,rep2.forceFilt(humerusStartIndex+i-1,:)); % N-mm resultant torque at the wrist given a force is applied at the ball and socket joint
%     %rep3.extraTorque = cross(ballPos'-wristPos,rep3.forceFilt(humerusStartIndex+i-1,:)); % N-mm resultant torque at the wrist given a force is applied at the ball and socket joint
%     rep1.saveExtraTorque(i,:) = rep1.extraTorque;
%     rep2.saveExtraTorque(i,:) = rep2.extraTorque;
%     %rep3.saveExtraTorque(i,:) = rep3.extraTorque;
%     rep1.torqueFilt(i,:) = (Jall(:,1:4)'*[rep1.forceFilt(humerusStartIndex+i-1,:)';rep1.extraTorque'])'/1000; % N-m
%     rep2.torqueFilt(i,:) = (Jall(:,1:4)'*[rep2.forceFilt(humerusStartIndex+i-1,:)';rep2.extraTorque'])'/1000; % N-m
%     %rep3.torqueFilt(i,:) = (Jall(:,1:4)'*[rep3.forceFilt(humerusStartIndex+i-1,:)';rep3.extraTorque'])'/1000; % N-m
% end
% if plots == 2
%     figure(10)
%     hold on
%     plot(keepTime,rep1.torqueFilt,'LineWidth',4)
%     plot(keepTime,rep2.torqueFilt,'LineWidth',2)
%     %plot(keepTime,rep3.torqueFilt,'LineWidth',1)
%     %plot(keepTime,rep1.torques(humerusStartIndex:humerusEndIndex,:),'o','LineWidth',4)
%     %plot(keepTime,rep2.torques(humerusStartIndex:humerusEndIndex,:),'o','LineWidth',2)
%     %plot(keepTime,rep3.torques(humerusStartIndex:humerusEndIndex,:),'o','LineWidth',1)
%     
%     figure(13)
%     hold on
%     plot(keepTime,rep1.forceFilt(humerusStartIndex:humerusEndIndex,:),'LineWidth',4)
%     plot(keepTime,rep2.forceFilt(humerusStartIndex:humerusEndIndex,:),'LineWidth',2)
%     %plot(keepTime,rep3.forceFilt(humerusStartIndex:humerusEndIndex,:),'LineWidth',1)
%   
% end
% 
% fill the initial times for which the humerus rigid body was not visable
% with the first set of joint angles for which is was
filteredJointAngles = [ones(humerusStartIndex-1,1)*[Xsmooth(10:12,1)' Xsmooth(14:16,1)'];
                       Xsmooth(10:12,:)' Xsmooth(14:16,:)';
                       ones(length(time)-humerusEndIndex,1)*[Xsmooth(10:12,end)' Xsmooth(14:16,end)']];
filteredWristPosition = [ones(humerusStartIndex-1,1)*wrist(1,:);
                        wrist;
                       ones(length(time)-humerusEndIndex,1)*wrist(end,:)];                   
% velocities      = [Xsmooth(map.aD(5):map.aD(7),:)' Xsmooth(map.aD(9):map.aD(11),:)'];
% accelerations   = [Xsmooth(map.aDD(5):map.aDD(7),:)' Xsmooth(map.aDD(9):map.aDD(11),:)'];


