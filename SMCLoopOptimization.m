%Script to optimize by having multiple SMC rounds, updating the prior at
%every step by adding the distribution from the last SMC

clear

%%
%Computation of Y_data

timeVector=0:0.5:10;%times at which measurements are taken

x0=[1,0,0,0,1,0];%Initial conditions. THIS CHANGES DEPENDING ON THE PROBLEM
kdeg=0.5;%degradation constant
rmax=2;%maximal rate of production
km=0.5;%Michaelis-Menten constant

%%%%THIS CHANGES DEPENDING ON THE PROBLEM
path1=4;
path2=2;
trueStrength1=0;
trueStrength2=0;
trueStrength3=1;
observableNode=6;
%%%%Simulation of the data generation (one trajectory)

Y=feval('twoPaths3',timeVector,x0,[kdeg,rmax,km,path1,path2,trueStrength1,trueStrength2,trueStrength3]);
Y_data=Y.statevalues(:,observableNode)+0.05*randn(length(Y.statevalues(:,observableNode)),1);%this variable contains what we actually can observe experimentally (additive gaussian noise added)

%%
meanTheta=[0,0,0;0.7,0.7,0.7];
totalMarginalLikelihood=[];
for j=1:5
%     meanTheta=[meanTheta(1,:);1,1,1];
    for i=1:5
        i
        [postTheta,marg_likelihood]=SMC_toyMEXadaptiveFixedParametersLoopOptimization(meanTheta,Y_data);
        meanTheta(1,:)=mean(log(postTheta([1,2,3],:)),2);
%         meanTheta(2,1)=sqrt(var(log(postTheta(1,:))));
%         meanTheta(2,2)=sqrt(var(log(postTheta(2,:))));
%         meanTheta(2,3)=sqrt(var(log(postTheta(3,:))));

    %     meanTheta(2,1)=sqrt(0.5*log(var(postTheta(1,:))));
    %     meanTheta(2,2)=sqrt(0.5*log(var(postTheta(2,:))));
    %     meanTheta(1,:)=2*log(mean(postTheta([1,2],:),2))./meanTheta(2,:).^2;
        totalMarginalLikelihood=[totalMarginalLikelihood;marg_likelihood];
    end
end
plot(log(totalMarginalLikelihood))