function A = mean5(hyp, x, i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEANMUSCLE.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Simple mean function for the toque about the shoulder that 
% results from stimulating a muscle. This function is called by the gp() 
% function of the GPML toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% hyp = [constantTerm thetaTerm];
%
% Inputs:
% x(:,1) = elevation plane 
% x(:,2) = shoulder elevation
% x(:,3) = shoulder rotation
% x(:,4) = elbow flexion
% x(:,5) = elbow pronation
%
% i = index of the hyperparameter which you want to take a derivative with
%     respect to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15 Febraury 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See also MEANFUNCTIONS.M.

if nargin<2, A = '2'; return; end             % report number of hyperparameters 

% parameter vector
params = hyp;
n = size(x,1);

PHI = compute_phi5(x)';

if nargin==2  % evaluate mean
    A = PHI*params;    % this is an n X 1 vector of means                                          
else % evaluate derivative (all zeros because there are no parameters 
    A = zeros(n,1); % this is an n X 1 vector of derivatives wrt the ith parameter 
end