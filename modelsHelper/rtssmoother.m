function [X, P] = rtssmoother(Xminus,Xplus,Pminus,Pplus,A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RTSSMOOTHER.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Descritpion: This function does RTS smoothing given that you have already
% done the forward filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   Xminus  = Dxn matrix of prior state estimates from the forward EKF
%   Xplus   = Dxn matrix of posterior state estimates from the forward EKF
%   Pminus  = DxDxn matrix of prior covariance estimates from the forward EKF
%   Pplus   = DxDxn matrix of posterior covariance estimates from the forward EKF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%   X       = Dxn matrix of smoothed state estimates 
%   P       = DxDxn matrix of smoothed covariance estimates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 30 July 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nD = size(Xplus,2);                                             % length of data sequence

X = zeros(size(Xminus,1),size(Xminus,2));                       % allocate space for state estimates
P = zeros(size(Pminus,1),size(Pminus,2),size(Pminus,3));        % allocate space for state estimates

% the final point in the sequence is just the forward posterior estimate 
X(:,end) = Xplus(:,end);
P(:,:,end) = Pplus(:,:,end);
for i = nD-1:-1:1
    scriptI = inv(Pminus(:,:,i+1));                                 % inverse of covariance
    K = Pplus(:,:,i)*A'*scriptI;                                    % smoother gain
    P(:,:,i) = Pplus(:,:,i) - K*(Pminus(:,:,i+1)-P(:,:,i+1))*K';    % covariance update
    X(:,i) = Xplus(:,i) + K*(X(:,i+1)-Xminus(:,i+1));               % state update
end