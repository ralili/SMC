%This script contains the code to perform the least-squares identification
%of the parameters (Theta). It depends on toyModelMinimize function.

numCycles=100;
res=zeros(numCycles,2);
for i=1:numCycles
    timeVector=0:0.5:10;%times at which measurements are taken
    x0=[1,0,0,1,0];%Initial conditions. THIS CHANGES DEPENDING ON THE PROBLEM
    kdeg=0.5;%degradation constant
    rmax=2;%maximal rate of production
    km=0.5;%Michaelis-Menten constant
    path1=3;
    path2=2;
    observableNode=5;
    trueStrength1=0;
    trueStrength2=1;

    Y=feval('twoPaths2',timeVector,x0,[kdeg,rmax,km,path1,path2,trueStrength1,trueStrength2]);
    Y_data=Y.statevalues(:,observableNode)+0.05*randn(length(Y.statevalues(:,observableNode)),1);

    res(i,:)=lsqnonlin(@(input) toyModelMinimize(input,Y_data),[0.5,0.5]);
    
end
mean(res,1)