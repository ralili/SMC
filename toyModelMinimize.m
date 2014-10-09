%This function is used to find the error between the input and the data. It
%is used to find the least-squares estimate of the parameters (Theta).

function [error]=toyModelMinimize(input,Y_data)

    timeVector=0:0.5:10;%times at which measurements are taken
    x0=[1,0,0,1,0];%Initial conditions. THIS CHANGES DEPENDING ON THE PROBLEM
    kdeg=0.5;%degradation constant
    rmax=2;%maximal rate of production
    km=0.5;%Michaelis-Menten constant
    path1=3;
    path2=2;
    observableNode=5;
    P_data=feval('twoPaths2',timeVector,x0,[kdeg,rmax,km,path1,path2,input(1),input(2)]);
    error=0.5*sum((Y_data-P_data.statevalues(:,observableNode)).^2);