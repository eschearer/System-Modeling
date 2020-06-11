function trainOut=TrainWristGP_multHyp(DoF,muscle,trainingInputs,trainingOutputs,params)
%%%%%%%%
hyp.cov=params;
hyp.lik=log(10);
%hyp.mean=[0 0];
hyp.mean=[]; % For no mean func

if strcmp(DoF,'Fx')
    meanfunc = @mean1;
end
if strcmp(DoF,'Fy')
    meanfunc = @mean2;
end
if strcmp(DoF,'Fz')
    meanfunc = @mean3;
end

meanfunc=[]; % for no mean func
inffunc={@infExact};
likfunc = {@likGauss}; 
covfunc={@covSE_rigid_multHyp};


hyp2 = minimize(hyp, @gp, -100, inffunc, meanfunc, covfunc, likfunc, trainingInputs, trainingOutputs);
[mtrain,s2train] = gp(hyp2, inffunc, meanfunc, covfunc, likfunc, trainingInputs, trainingOutputs, trainingInputs);

% training error and RMS
trainingError1 = trainingOutputs(:,1)-mtrain;
trainingRMS1 = sqrt(mean(trainingError1.^2));

h1=figure(1);
plot(mtrain,trainingOutputs,'*')
saveas(h1,['./figures/',DoF,'Muscle',num2str(muscle),'TestVsTrain.jpg'])
saveas(h1,['./figures/',DoF,'Muscle',num2str(muscle),'TestVsTrain.fig'])
drawnow
h2 = figure(2);
z = linspace(1,length(mtrain),length(mtrain));
hv = fill([z'; flipdim(z',1)],[mtrain+2*sqrt(s2train);flipdim(mtrain-2*sqrt(s2train),1)],[7 7 7]/8);
hold on
htSIMPLET1 =plot(trainingOutputs(:,1),'b','LineWidth',2);
hold off
ylabel('Force')
saveas(h2,['./figures/',DoF,'Muscle',num2str(muscle),'nonParametric.jpg'])
saveas(h2,['./figures/',DoF,'Muscle',num2str(muscle),'nonParametric.fig'])




hyp2


trainOut.trainingError        = trainingError1;
trainOut.trainingRMS          = trainingRMS1;
trainOut.trainingPredictions  = mtrain;
trainOut.trainingInputs       = trainingInputs;
trainOut.trainingOutputs      = trainingOutputs;
trainOut.hyp                  = hyp2;
end