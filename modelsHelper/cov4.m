function K = cov4(hyp, x, z, i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COVMUSCLE1.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Simple cov function for the torque about the shoulder that 
% results from stimulating a muscle. This function is called by the gp() 
% function of the GPML toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hyp = [  
%          B ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15 February 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See also COVFUNCTIONS.

if nargin<2, K = '3'; return; end                  % report number of hyperparameters 

% this is the number of training samples
[n,D] = size(x);

% here are the basis functions
PHIX = compute_phi4(x)';

if nargin<3, z = []; end                                   % make sure z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

% this is the covariance matrix of the parameters
B = [hyp(1) hyp(2);  
     hyp(2) hyp(3)]; 

%%%%%%%%%%%%%%%%%%%%%%%%
phiPart = zeros(n,1);
if nargin<4                                                        % covariances
  if dg
      for j = 1:n
          phiPart(j) =  PHIX(j,:)*B*PHIX(j,:)';
      end
      K = phiPart; % vector kxx
  else
      if xeqz
          K = PHIX*B*PHIX';  % symmetric matrix Kxx
      else
          nz = size(z,1);
          PHIZ = compute_phi4(z)';
          K = PHIX*B*PHIZ';  % cross covariances Kxz
      end
  end 
else          % derivatives
    if i>0 && i <= 21
        if dg
            K = zeros(n,1);
        else
            if xeqz
                K = zeros(n);
            else
                K = zeros(n,nz);
            end
        end
    else
        error('Unknown hyperparameter')
    end
end

