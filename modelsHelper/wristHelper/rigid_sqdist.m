%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rigid_sqdist(q1,q2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the squared distance between 2 rigid bodies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%       q1: matrix of row vectors of position and quaternion for body 1 nx7
%       q2: matrix of row vectors of position and quaternion for body 2 mx7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Derek Wolf   3/21/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revisions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sqdist] = rigid_sqdist(q1,q2)
% Weighting
rho1=1;
rho2=.1;

% Select data
pos1=q1(:,1:3);
pos2=q2(:,1:3);
quat1=q1(:,4:7);
quat2=q2(:,4:7);

% Get axis angle
[~,axis1,angle1] = quat2axisAngle2(quat1);
[~,axis2,angle2] = quat2axisAngle2(quat2);

n=size(q1,1);
m=size(q2,1);


% Arc distance
for i=1:n
    for j=1:m
        d_arc(i,j)=2.*acos(abs(cos(angle1(i)/2)*cos(angle2(j)/2)+sin(angle1(i)/2)*sin(angle2(j)/2)*axis1(i,:)*axis2(j,:)'));
    end
end
d_arc=d_arc.^2;

% Vector distance
v_dist=sq_dist(pos1',pos2');

sqdist=rho1.*v_dist+rho2.*d_arc;