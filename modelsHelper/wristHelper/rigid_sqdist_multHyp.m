%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rigid_sqdist(q1,q2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the squared distance between 2 rigid bodies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%       q1: matrix of row vectors of position and quaternion for body 1 nx7
%       q2: matrix of row vectors of position and quaternion for body 2 mx7
%       lam: vector of lambdas
%           1 2 3 are the vector distance
%           4 is the arc distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Derek Wolf   3/21/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revisions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sqdist,d_arc,v_dist,v_dist_der] = rigid_sqdist_multHyp(q1,q2,lam,der_lam)
% Weighting
rho1=1; %.1; %1;
rho2=0.1; %0.1%

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
d_arc=zeros(n,m);
for i=1:n
    for j=1:m
        d_arc(i,j)=2.*acos(abs(cos(angle1(i)/2)*cos(angle2(j)/2)+sin(angle1(i)/2)*sin(angle2(j)/2)*axis1(i,:)*axis2(j,:)'));
    end
end
d_arc=rho2.*d_arc.^2/(2*lam(4)^2);

% Vector distance
lamdiag=diag([1/lam(1),1/lam(2),1/lam(3)]);
lamdiag=lamdiag.^2/2;
v_dist=zeros(n,m);
for i=1:n
    for j=1:m
        v_dist(i,j)=(pos1(i,:)'-pos2(j,:)')'*lamdiag*(pos1(i,:)'-pos2(j,:)');
    end
end
%v_dist=sq_dist(pos1',pos2');
v_dist=rho1*v_dist;

sqdist=v_dist+d_arc;


% v_dist_der
v_dist_der=zeros(n,m);
if nargin==4
    for i=1:n
        for j=1:m
            v_dist_der(i,j)=rho1*(pos1(i,1)-pos2(j,1))^2/(lam(der_lam)^3);
        end
    end
end














