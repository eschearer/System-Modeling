function K = covSE_rigid_multHyp(hyp, x, z, i)

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



if nargin<2, K = '5'; return; end              % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode


sf=hyp(1);
lam1=hyp(2);
lam2=hyp(3);
lam3=hyp(4);
lam4=hyp(5); % d_arc

lam=[lam1 lam2 lam3 lam4];


if dg
    K=zeros(size(x,1),1);
else
    if xeqz
        K=rigid_sqdist_multHyp(x,x,lam);
    else
        
        K=rigid_sqdist_multHyp(x,z,lam);
    end
end
K = sf^2*exp(-K);

if nargin>3
    if i==1 % derivative wrt sf
        if dg
            K=K*0;
            size(K)
            pause
        else
            K=2*K/sf;
        end
    elseif i==2 % derivative wrt lam1
        lamID=i-1;
        if xeqz
            [~,~,~,der]=rigid_sqdist_multHyp(x,x,lam,lamID);
        else
            [~,~,~,der]=rigid_sqdist_multHyp(x,z,lam,lamID);
        end
        K=der.*K;
    elseif i==3 % derivative wrt lam2
        lamID=i-1;
        if xeqz
            [~,~,~,der]=rigid_sqdist_multHyp(x,x,lam,lamID);
        else
            [~,~,~,der]=rigid_sqdist_multHyp(x,z,lam,lamID);
        end
        K=der.*K;
    elseif i==4 % derivative wrt lam3
        lamID=i-1;
        if xeqz
            [~,~,~,der]=rigid_sqdist_multHyp(x,x,lam,lamID);
        else
            [~,~,~,der]=rigid_sqdist_multHyp(x,z,lam,lamID);
        end
        K=der.*K;
        
    elseif i==5 % derivative wrt lam4
        if xeqz
            [~,der,~,~]=rigid_sqdist_multHyp(x,x,lam);
        else
            [~,der,~,~]=rigid_sqdist_multHyp(x,z,lam);
        end
        der=der.*2/lam4; % = rho2*d_arc^2/lam(4)^3
        K=der.*K;
    else
        error('Unknown hyperparameter');
    end
end
end

