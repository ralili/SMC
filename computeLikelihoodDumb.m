result4=zeros(40,1);
for i=1:40
    meanlogpar=[-2,-2,-2,0,-2,-2,-2];

    loglikelihood = zeros(blength,M);  % likelihood before resampling
    sigmanoise = 0.05;

    sigmalogpar=ones(1,length(meanlogpar));
    timeVector=0:0.5:10;%times at which measurements are taken

    x0=[1,0,0,0,1,0,0,0];%Initial conditions
    kdeg=0.5;%degradation constant
    kmax=2;%maximal rate of production
    km=0.5;%Michaelis-Menten constant
    connectivity=[1 5;1 6;1 8;2 5;2 6;2 7;4 8];
    trueStrength=[1;0;0;1;0;0;0];
    [~, Y] = ode45(@(t,x) toyModel(t,x,kdeg,kmax,km,connectivity,trueStrength),timeVector,x0);%Simulation of real system
    Y_data=Y(:,8)+0.05*randn(length(Y(:,8)),1);%this variable contains what we actually can observe experimentally (additive gaussian noise added)

    parnum=length(trueStrength);%number of parameters

    initial_points = zeros(parnum,M);
    for k = 1:parnum
        initial_points(k,:) = lognrnd(meanlogpar(k),sigmalogpar(k),M,1);  % draw initial points distributed according to prior (log-normal)
    end

    for m = 1:M
        [~, P_data] = ode45(@(t,x) toyModel(t,x,kdeg,kmax,km,connectivity,initial_points(:,m)),timeVector,x0);  %simulating data from each particle
        loglikelihood(1,m) = 0.5*sum((Y_data-P_data(:,8)).^2)/sigmanoise;                   %likelihood computation
    end

    result4(i)=mean(loglikelihood(1,:));
    if mean(loglikelihood(1,:))<150
        mean(initial_points(:,:),2)
    end
 end