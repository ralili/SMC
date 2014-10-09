%

clear
clf

M = 200;
sigmanoise=0.05;

timeVector=0:0.5:10;%times at which measurements are taken

x0=[1,0,0,1,0];%Initial conditions. THIS CHANGES DEPENDING ON THE PROBLEM
kdeg=0.5;%degradation constant
rmax=2;%maximal rate of production
km=0.5;%Michaelis-Menten constant

%%Here the known model parameters are given prior distributions
variable_parameters=[1,2];
%%%THIS CHANGES DEPENDING ON THE PROBLEM
path1=3;
path2=2;
trueStrength1=0;
trueStrength2=1;
observableNode=5;

Y=feval('twoPaths2',timeVector,x0,[km,rmax,kdeg,path1,path2,trueStrength1,trueStrength2]);
Y_data=Y.statevalues(:,observableNode)+0.05*randn(length(Y.statevalues(:,observableNode)),1);%this variable contains what we actually can observe experimentally (additive gaussian noise added)


%

hyperparameters=-4.8:0.4:1;
likelihood=zeros(length(hyperparameters));

for i=1:length(hyperparameters)
    i
    for j=1:length(hyperparameters)
        loglikelihood=zeros(1,M);
        meanlogpar=[hyperparameters(i),hyperparameters(j)];
        sigmalogpar=ones(1,length(meanlogpar)).*[0.7,0.01];
%         sigmalogpar=0.1./exp(meanlogpar);
        meanlogpar=[meanlogpar,log([km,rmax,kdeg])];
        sigmalogpar=[sigmalogpar,0,0,0];%%%%%%
        parnum=length(meanlogpar);

        parfor k = 1:parnum
            initial_points(k,:) = lognrnd(meanlogpar(k),sigmalogpar(k),M,1);  % draw initial points distributed according to prior (log-normal)
        end

        parfor m = 1:M
            P_data=feval('twoPaths2',timeVector,x0,[initial_points(end-2:end,m)',path1,path2,initial_points(1:2,m)']);
            loglikelihood(m) = 0.5*sum((Y_data-P_data.statevalues(:,observableNode)).^2)/sigmanoise;                   %likelihood computation
        end
        likelihood(i,j)=mean(loglikelihood);
    end
end

likelihood2=zeros(length(hyperparameters));%Likelihood computed just looking at a certain parameter value, and not using a prior distribution.

for i=1:length(hyperparameters)
    for j=1:length(hyperparameters)
        P_data=feval('twoPaths2',timeVector,x0,[km,rmax,kdeg,path1,path2,exp(hyperparameters(i)),exp(hyperparameters(j))]);
        likelihood2(i,j) = 0.5*sum((Y_data-P_data.statevalues(:,observableNode)).^2)/sigmanoise;                   %likelihood computation
    end
end

figure(1)
subplot(1,2,1)
imagesc(hyperparameters,hyperparameters,likelihood,[min(likelihood(:)),1e2])
xlabel('Strength2, logarithmic')
ylabel('Strength1li, logarithmic')
colorbar
set(gca,'YDir','normal')
title('Hyper')
subplot(1,2,2)
imagesc(hyperparameters,hyperparameters,likelihood2,[min(likelihood2(:)),1e2])
colorbar
set(gca,'YDir','normal')
xlabel('Strength2, logarithmic')
ylabel('Strength1, logarithmic')
title('Likelihood')

figure(2)
subplot(1,2,1)
[minimum,index]=min(likelihood(:));

[I,J] = ind2sub(size(likelihood),index);
OptStrength1=hyperparameters(I);
OptStrength2=hyperparameters(J);

P_data=feval('twoPaths2',timeVector,x0,[0.5,2,0.5,path1,path2,exp(OptStrength1),exp(OptStrength2)]);

hold all
plot(Y_data)
plot(P_data.statevalues(:,observableNode))

disp(OptStrength1)
disp(OptStrength2)
hold off

subplot(1,2,2)

[minimum,index]=min(likelihood2(:));

[I,J] = ind2sub(size(likelihood2),index);
OptStrength1=hyperparameters(I);
OptStrength2=hyperparameters(J);

P_data=feval('twoPaths2',timeVector,x0,[0.5,2,0.5,path1,path2,exp(OptStrength1),exp(OptStrength2)]);

hold all
plot(Y_data)
plot(P_data.statevalues(:,observableNode))


disp(OptStrength1)
disp(OptStrength2)
hold off
