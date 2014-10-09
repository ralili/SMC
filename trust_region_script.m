%This script is used to optimize the inputs for a function, provided that
%the function ouputs the Gradient and the Hessian of the objective function
%(likelihood).

clear

%%
%Computation of Y_data

timeVector=0:0.5:10;%times at which measurements are taken

x0=[1,0,0,1,0];%Initial conditions. THIS CHANGES DEPENDING ON THE PROBLEM
kdeg=0.5;%degradation constant
rmax=2;%maximal rate of production
km=0.5;%Michaelis-Menten constant

%%%%THIS CHANGES DEPENDING ON THE PROBLEM
path1=3;
path2=2;
trueStrength1=0;
trueStrength2=1;
observableNode=5;
%%%%Simulation of the data generation (one trajectory)

Y=feval('twoPaths2',timeVector,x0,[kdeg,rmax,km,path1,path2,trueStrength1,trueStrength2]);
Y_data=Y.statevalues(:,observableNode)+0.05*randn(length(Y.statevalues(:,observableNode)),1);%this variable contains what we actually can observe experimentally (additive gaussian noise added)
%%
%Optimization procedure

f = @(hyperparameters) SMC_toyMEXadaptiveFixedParameters(hyperparameters,Y_data);
hyperparameters0=[0,-6;0,0];

options = optimset('GradObj','on','Hessian','on','Display','iter','MaxFunEvals',4,'TolFun',1e-7);
[meanlogpar,fval,exitflag,output] = fminunc(f,hyperparameters0,options);

%SMC_toyMEXadaptiveFixedParameters(meanlogpar,Y_data);
