function K = covSE_rigid(hyp, x, z, i)

% Squared exponential covariance function for rigid body transformations
% K = sf2 * exp(-dist^2/(2*lam))
% parameters are sf2 and lam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Derek Wolf 3/21/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revisions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Built based on the following toolbox
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.



if nargin<2, K = '2'; return; end              % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode



[~,D] = size(x);

sf=hyp(1);
lam=hyp(2);




if dg
    K=zeros(size(x,1),1);
else
    if xeqz
        K=rigid_sqdist(x,x);
    else
        
        K=rigid_sqdist(x,z);
    end
end
K = sf^2*exp(-K./(2*lam^2));

if nargin>3
    if i==1
        if dg
            K=K*0;
            size(K)
            pause
        else
            K=2*K/sf;
            
        end
    elseif i==2
        if xeqz
            der=rigid_sqdist(x,x);
        else
            
            der=rigid_sqdist(x,z);
        end
        der=der/(lam^3);
        K=K.*der;
    else
        error('Unknown hyperparameter');
    end
end
