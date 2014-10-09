%%
%Generation of experimental data

x0=[1,0,0,0,1,0,0,0];%Initial conditions
kdeg=0.5;%degradation constant
kmax=2;%maximal rate of production
km=0.5;%Michaelis-Menten constant
[T, Y] = ode45(@(t,x) toyModel(t,x,kdeg,kmax,km,[2 6;2 7;2 8],[0 1 0]),0:10,x0);%Simulation of real system
observed=Y(:,8)+0.05*randn(length(Y(:,8)),1);%this variable contains what we actually can observe experimentally (additive gaussian noise added)

%%
%Example of how to compute the likelihood for various models

[T1, model1] = ode45(@(t,x) toyModel(t,x,kdeg,kmax,km,[2 6;2 7;2 8],[1 0 0]),0:10,x0);%Simulation
likelihood1=sum(-(observed-model1(:,8))'*(0.05)^-1*(observed-model1(:,8)));

[T2, model2] = ode45(@(t,x) toyModel(t,x,kdeg,kmax,km,[2 6;2 7;2 8],[0 1 0]),0:10,x0);%Simulation
likelihood2=sum(-(observed-model2(:,8))'*(0.05)^-1*(observed-model2(:,8)));

[T3, model3] = ode45(@(t,x) toyModel(t,x,kdeg,kmax,km,[2 6;2 7;2 8],[0 0 1]),0:10,x0);%Simulation
likelihood3=sum(-(observed-model3(:,8))'*(0.05)^-1*(observed-model3(:,8)));

%%
%Example of how to marginalize over theta, given certain priors, to obtain
%the marginalized likelihood. Very rudamentary implementation, using Monte
%Carlo integration.

likelihood=0;
for i=1:1000
    theta=[exp(-10000+randn(1)) exp(0+randn(1)) exp(-10000+randn(1))];
    [~,temp] = ode45(@(t,x) toyModel(t,x,kdeg,kmax,km,[2 6;2 7;2 8],theta),0:10,x0);
    likelihood=likelihood+sum(-(observed-temp(:,8))'*(0.05)^-1*(observed-temp(:,8)));
end
integral=likelihood/1000;

%%
%Computing P(D|prior) via sequential Monte Carlo with Annealed Importance
%Sampling.

%Sample the parameter vector from prior distribution and compute the
%likelihood of those samples.

%Compute the distribution of the resulting points via a Gaussian mixture

%Perform a few steps of MCMC with the Gaussian mixture as proposal
%distribution and prior^beta*P(theta|D)^(1-beta) as target distribution.
%*How do we know what P(theta|D) looks like? Likelihood function? Not
%really

%Set adaptive beta

