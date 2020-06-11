%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTEFEASIBLE_wrist.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: This script checks to see which static wrist positions can
% be achieved given the model of a particular subject on a particular day.
% It outputs a file called 'feasiblepoints.mat' that contains the feasible
% wrist positions, the associated joint angles, and the muscle activations
% required to hold the positions.  This is called "Enhanced" because it
% allows for increasing the strength of an individual muscle.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 24 July 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated:
%   15 Feb 2016:    Made updates to use the new data format for November
%                   2015 experiments.
%   9 Feb 2017:     Added a flag to downsample the data to make GP run
%                   faster
%   7 Sep 2017:     Added inverse kinematics, minimizing the required
%                   torque and equal grid spacing in cartesian space
%                   instead of joint space.
%   3/22/18:        Made it based on wrist positions not joint angle
%                   control - Derek Wolf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function computeFeasible_wrist(filename)

close all

global modeldata

sampling = 1; % only use 1/5 of the data to fit models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load muscle models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set these constants based on initial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a set of joint angles based on the min and max joint angles from the
% identification experiments and whether those joint angles are inside the
% convex hull of actual wrist positions visited
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visited wrist positions in Thorax Coordinate Frame
visitedWristPose =   [unique(modeldata.muscle(1).Fx.trainingInputs,'rows');
    unique(modeldata.muscle(2).Fx.trainingInputs,'rows');
    unique(modeldata.muscle(3).Fx.trainingInputs,'rows');
    unique(modeldata.muscle(4).Fx.trainingInputs,'rows');
    unique(modeldata.muscle(5).Fx.trainingInputs,'rows');
    unique(modeldata.muscle(6).Fx.trainingInputs,'rows');
    unique(modeldata.muscle(7).Fx.trainingInputs,'rows');
    unique(modeldata.muscle(8).Fx.trainingInputs,'rows');
    unique(modeldata.muscle(9).Fx.trainingInputs,'rows');
    unique(modeldata.muscle(10).Fx.trainingInputs,'rows')];

%% create a grid of positions inside the convex hull of visited wrist positions
visitedWristPositions=visitedWristPose(:,1:3);
minPositions = min(visitedWristPositions);
maxPositions = max(visitedWristPositions);
spacing = 10; % grid spacing in mm
xCoords = minPositions(1):spacing:maxPositions(1);
yCoords = minPositions(2):spacing:maxPositions(2);
zCoords = minPositions(3):spacing:maxPositions(3);
[X,Y,Z] = meshgrid(xCoords,yCoords,zCoords);
X = reshape(X,length(xCoords)*length(yCoords)*length(zCoords),1);
Y = reshape(Y,length(xCoords)*length(yCoords)*length(zCoords),1);
Z = reshape(Z,length(xCoords)*length(yCoords)*length(zCoords),1);
inhullIndices = find(inhull([X Y Z],visitedWristPositions));
testWristPositions = [X(inhullIndices) Y(inhullIndices) Z(inhullIndices)];
plot(testWristPositions(:,1),-testWristPositions(:,3),'o')
hold on
plot(visitedWristPositions(:,1),-visitedWristPositions(:,3),'x')
axis image

%%
%%%%%%%%%%%%%%%%%%
% Find nearest neighbor of visited wristPosition and use for orientation
idx=knnsearch(visitedWristPositions,testWristPositions);
visitedWristQuaternion=visitedWristPose(:,4:7);

% create test vector
q=[testWristPositions visitedWristQuaternion(idx,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comment out to do actual test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%q = modeldata.muscle(10).Fx.trainingInputs(1:sampling:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through the joint angles and determine the stimulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%t = zeros(size(q,1),3);
activation = zeros(size(q,1),9);
hmForce = zeros(size(q,1),3);
stimulation = zeros(size(q,1),9);
flag = zeros(size(q,1),1);

meanfunc=[];
likfunc = {@likGauss};
covfunc={@covSE_rigid_multHyp};
%covfunc={@covSE_rigid};

% modeldata.muscle(10).Fx.hyp.cov(2)=20;
% modeldata.muscle(10).Fy.hyp.cov(2)=20;
% modeldata.muscle(10).Fz.hyp.cov(2)=20;
t(:,1) = gp(modeldata.muscle(10).Fx.hyp, @infExact, meanfunc, covfunc, likfunc, modeldata.muscle(10).Fx.trainingInputs(1:sampling:end,:), modeldata.muscle(10).Fx.trainingOutputs(1:sampling:end), q);
t(:,2) = gp(modeldata.muscle(10).Fy.hyp, @infExact, meanfunc, covfunc, likfunc, modeldata.muscle(10).Fy.trainingInputs(1:sampling:end,:), modeldata.muscle(10).Fy.trainingOutputs(1:sampling:end), q);
t(:,3) = gp(modeldata.muscle(10).Fz.hyp, @infExact, meanfunc, covfunc, likfunc, modeldata.muscle(10).Fz.trainingInputs(1:sampling:end,:), modeldata.muscle(10).Fz.trainingOutputs(1:sampling:end), q);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the passive equilibrium configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[val,startIndex] = min(abs(t(:,1)) + abs(t(:,2)) + abs(t(:,3)));
fun = @requiredforce;
options = optimset('Display','iter');
qEq = fmincon(fun,q(startIndex,:),[],[],[],[],min(q),max(q),[]);

% add the passive equilibrium configuration to q and make the corresponding
% torque zero
q(length(q)+1,:) = qEq;
t(length(q),:) = [0 0 0];
display('found equilibrium configuration')
%%
M = zeros(3,9,length(q));
sequence = [1 2 3 4 5 6 7 8 9];


for k = 1:length(sequence)
    k
    j = sequence(k);
    % Selected by trial and error
%     modeldata.muscle(j).Fx.hyp.cov(2)=20;
%     modeldata.muscle(j).Fy.hyp.cov(2)=20;
%     modeldata.muscle(j).Fz.hyp.cov(2)=20;
    
    
    M(1,j,:) = t(:,1) - gp(modeldata.muscle(j).Fx.hyp, @infExact, meanfunc, covfunc, likfunc, modeldata.muscle(j).Fx.trainingInputs(1:sampling:end,:), modeldata.muscle(j).Fx.trainingOutputs(1:sampling:end), q);
    M(2,j,:) = t(:,2) - gp(modeldata.muscle(j).Fy.hyp, @infExact, meanfunc, covfunc, likfunc, modeldata.muscle(j).Fy.trainingInputs(1:sampling:end,:), modeldata.muscle(j).Fy.trainingOutputs(1:sampling:end), q);
    M(3,j,:) = t(:,3) - gp(modeldata.muscle(j).Fz.hyp, @infExact, meanfunc, covfunc, likfunc, modeldata.muscle(j).Fz.trainingInputs(1:sampling:end,:), modeldata.muscle(j).Fz.trainingOutputs(1:sampling:end), q);
    
end
%%
feasible = zeros(length(q),1);
wristPosition = zeros(length(q),3);
c = zeros(length(q),3);
s = zeros(length(q),1);
for i = 1:length(q)
    display([num2str(i),' of ',num2str(length(q))]);
    
    % solve the optimization problem to find activations
    H = eye(9);
    f = zeros(9,1);
    [a val flag(i)] = quadprog(H,f,[],[],[squeeze(M(:,:,i))],[t(i,:)]',-0.0*ones(9,1),1*ones(9,1));
    if flag(i) == 1
        activation(i,:) = a';
        feasible(i) = 1;
    else
        activation(i,:) = zeros(1,9);
    end
    %     if flag(i) == -2
    %         activation(i,:) = zeros(1,9);
    %     else
    %         activation(i,:) = a';
    %         feasible(i) = 1;
    %     end
    
    % compute the wrist position and assign each wrist position a color
    % based on if it is feasible
    wristPosition(i,:) = q(i,1:3);
    if feasible(i)
        c(i,:) = [0 255 0];
        s(i) = 10;
    else
        c(i,:) = [255 0 0];
        s(i) = 1;
    end
    
end

figure
h=scatter(wristPosition(:,1),-wristPosition(:,3),s,c);
xlabel('x')
ylabel('-z')
figure
h=scatter(wristPosition(:,1),wristPosition(:,2),s,c);
xlabel('x')
ylabel('y')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('feasible.mat','q','wristPosition','feasible','activation','t','M')





