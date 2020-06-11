%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REQUIREDTORQUE.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Descritpion: This function computes the torque required to hold a desired
% arm configuration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 18 Feb 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function torqueMag = requiredtorque(q)

global modeldata

t = zeros(1,4);

covGP = {@covSEardNoDer};

covSIMPLE1 = {@cov1};
covfunc1 = {'covSum',{covGP,covSIMPLE1}};
t(1) = gp(modeldata.muscle(10).elevationplane.semiparametric.gpHyperparameters, @infExact, @mean1, covfunc1, @likGauss, modeldata.muscle(10).elevationplane.semiparametric.trainingInputs, modeldata.muscle(10).elevationplane.semiparametric.trainingOutputs, q);

covSIMPLE2 = {@cov2};
covfunc2 = {'covSum',{covGP,covSIMPLE2}};
t(2) = gp(modeldata.muscle(10).shoulderelevation.semiparametric.gpHyperparameters, @infExact, @mean2, covfunc2, @likGauss, modeldata.muscle(10).shoulderelevation.semiparametric.trainingInputs, modeldata.muscle(10).shoulderelevation.semiparametric.trainingOutputs, q);

covSIMPLE3 = {@cov3};
covfunc3 = {'covSum',{covGP,covSIMPLE3}};
t(3) = gp(modeldata.muscle(10).shoulderrotation.semiparametric.gpHyperparameters, @infExact, @mean3, covfunc3, @likGauss, modeldata.muscle(10).shoulderrotation.semiparametric.trainingInputs, modeldata.muscle(10).shoulderrotation.semiparametric.trainingOutputs, q);

covSIMPLE4 = {@cov4};
covfunc4 = {'covSum',{covGP,covSIMPLE4}}; 
t(4) = gp(modeldata.muscle(10).elbowflexion.semiparametric.gpHyperparameters, @infExact, @mean4, covfunc4, @likGauss, modeldata.muscle(10).elbowflexion.semiparametric.trainingInputs, modeldata.muscle(10).elbowflexion.semiparametric.trainingOutputs, q);

torqueMag = norm(t);