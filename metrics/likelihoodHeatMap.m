%%
%Parameter setup
timeVector=0:0.5:10;%times at which measurements are taken
M=50;%number of samples per each point in the heatMap
loglikelihood=zeros(1,M);%empty vector containing the loglikelihood for each sample
sigmanoise=0.05;%measurement noise

x0=[1,0,0,1,0];%Initial conditions. THIS CHANGES DEPENDING ON THE PROBLEM
kdeg=0.5;%degradation constant
kmax=2;%maximal rate of production
km=0.5;%Michaelis-Menten constant

pathLengths=[3,2];
connectivity=[2 5;3 5];%THIS CHANGES DEPENDING ON THE PROBLEM
trueStrength=[0;1];%THIS CHANGES DEPENDING ON THE PROBLEM
observableNode=5;%THIS CHANGES DEPENDING ON THE PROBLEM

%%
%Data generation
[~, Y] = ode45(@(t,x) toyModel(t,x,kdeg,kmax,km,pathLengths,connectivity,[1,1,0],[1,0,-1],[0,1,0],[0,5,10]),timeVector,[1,0,0,1,0]);%Simulation of real system
Y_data=Y(:,observableNode)+0.05*randn(length(Y(:,observableNode)),1);%this variable contains what we actually can observe experimentally (additive gaussian noise added)

%%

meanlogpar=[0,0,kdeg,kmax,km];
sigmalogpar=[0.2,0.2,0.2,0.2,0.1];
parnum=length(meanlogpar);

index=-3.5:0.5:1.5;
results=zeros(length(index));
for i=1:length(index)
    for j=1:length(index)
        meanlogpar=[index(i),index(j),kdeg,kmax,km];
        initial_points = zeros(parnum,M);
        for k = 1:parnum
            initial_points(k,:) = lognrnd(meanlogpar(k),sigmalogpar(k),M,1);  % draw initial points distributed according to prior (log-normal)
        end
        for m = 1:M
            [~, P_data] = ode45(@(t,x) toyModel(t,x,initial_points(end-2,m),initial_points(end-1,m),initial_points(end,m),pathLengths,connectivity,initial_points(1:end-3,m),[1,0,-1],[0,1,0],[0,5,10]),timeVector,[1,0,0,1,0]);
            loglikelihood(1,m) = 0.5*sum((Y_data-P_data(:,observableNode)).^2)/sigmanoise;                   %likelihood computation
        end
        results(i,j)=mean(loglikelihood(1,:));
    end
end