function [ev,gF,H] = SMC_toy(meanlogpar)%evidence, gradient of the evidence, hessian of the evidence
%% Setup sampler
M = 200;   % number of particles to be propagated
beta_schedule = ((1:30)/30).^3;%%it was ^3
blength = length(beta_schedule);
burn = 7;      % number of burn-in samples  (could vary with the value of beta)
Theta = cell(1,blength);  % particles from all generations
loglikelihood = zeros(blength,M);  % likelihood before resampling
W = zeros(blength,M);
w = zeros(blength,M);
W(1,:) = 1/M;
incr_w = zeros(1,M);
r_times = [];
prior = zeros(1,M);
sigmanoise = 0.05;
meanlogpar
sigmalogpar=ones(1,length(meanlogpar));

%% Parameter setup
timeVector=0:0.5:10;%times at which measurements are taken

x0=[1,0,0,1,0];%Initial conditions. THIS CHANGES DEPENDING ON THE PROBLEM
kdeg=0.5;%degradation constant
kmax=2;%maximal rate of production
km=0.5;%Michaelis-Menten constant

%%%Here the known model parameters are given prior distributions
meanlogpar=[meanlogpar,kdeg,kmax,km];
sigmalogpar=[sigmalogpar,0.5,0.5,0.5];
%%%
pathLengths=[3,2];
connectivity=[2 5;3 5];%THIS CHANGES DEPENDING ON THE PROBLEM
trueStrength=[0;1];%THIS CHANGES DEPENDING ON THE PROBLEM
observableNode=5;%THIS CHANGES DEPENDING ON THE PROBLEM
if length(x0) ~= sum(pathLengths)
    disp('error. Check number of species (initial conditions and pathLengths)')
    return
elseif length(connectivity) ~= length(trueStrength)
    disp('error. Check trueStrength')
    return
elseif length(connectivity) ~= length(meanlogpar)-3
    disp('error. Check input (meanlogpar)')
    return
end


%% Data generation
[~, Y] = ode45(@(t,x) toyModel(t,x,kdeg,kmax,km,pathLengths,connectivity,trueStrength,[1,1],[1,1],[0,timeVector(end)]),timeVector,x0);%Simulation of real system
Y_data=Y(:,observableNode)+0.05*randn(length(Y(:,observableNode)),1);%this variable contains what we actually can observe experimentally (additive gaussian noise added)

parnum=length(meanlogpar);%number of parameters

initial_points = zeros(parnum,M);
for k = 1:parnum
    initial_points(k,:) = lognrnd(meanlogpar(k),sigmalogpar(k),M,1);  % draw initial points distributed according to prior (log-normal)
end

%% GM setup
options = statset('MaxIter',500,'TolFun',1e-5);

%% Sampler

for m = 1:M
    [~, P_data] = ode45(@(t,x) toyModel(t,x,initial_points(end-2,m),initial_points(end-1,m),initial_points(end,m),pathLengths,connectivity,initial_points(1:end-3,m),[0,0],[0,0],[0,timeVector(end)]),timeVector,x0);  %simulating data from each particle
    loglikelihood(1,m) = 0.5*sum((Y_data-P_data(:,observableNode)).^2)/sigmanoise;                   %likelihood computation
end

w(2,:) = exp(-beta_schedule(1)*(loglikelihood(1,:)));                               %calculate incremental weights w_1(X_1^i,X_2^i) from {X_1^(i)}. Isn't the likelihood from the prior missing?
if (M/(1+var(w(2,:)/mean(w(2,:))))<0.5*M)           %If effect sample size is below 0.5, then resample. Resampling decreases weight variance, therefore not letting effect size go too low.
    resample = 1;
    r_times = 1;
    prob = w(2,:)/sum(w(2,:));
    W(2,:) = 1/M;                                                                      %calculate {W_2^i} (different if resampling occurs)
else
    resample = 0;
    W(2,:) = w(2,:)/sum(w(2,:));
    prob = 1/M;   %(just to define pr later)
end

Theta{1} = initial_points;                                                             %first generation {X_1^i} - before resampling


dispstat('','init');
dispstat('Running SMC...','timestamp','keepthis');


for b = 2:blength   % b = n (from paper)
    dispstat(sprintf('Step %d of 30',b),'timestamp');
    beta = beta_schedule(b-1);
    prev_theta = Theta{b-1};
    next_theta = zeros(parnum,M);
    if (resample)
        prob(find(prob==0)) = min(prob(find(prob~=0)))/1000;% What does this do?? Is it to pick up samples from the weighted particles?
        replicates = floor(prob*M);
        repl_ind = [];
        theta_res = [];
        for u = 1:M
            theta_res = [theta_res Theta{b-1}(:,u*ones(1,replicates(u)))];
            repl_ind = [repl_ind u*ones(1,replicates(u))];
        end
        rest = sum(replicates);
        for u = rest+1:M
            rndpoint = randsample(M,1,true,prob);
            theta_res = [theta_res Theta{b-1}(:,rndpoint)];
            repl_ind = [repl_ind rndpoint];
        end
        prev_theta = theta_res;
    end
    pr = prob;%What does this do?
    LL = loglikelihood(b-1,:);
    accepts = 0;
    clusters = 6;
    [GMMtheta,NlogL,optimInfo] = clusterGM(log(prev_theta)',clusters,'rndsample',1,2,0,1e-6,options);
    Chol_Sigma = zeros(size(GMMtheta.Sigma));
    for i = 1:clusters
        Chol_Sigma(:,:,i) = cholcov_A(GMMtheta.Sigma(:,:,i));   %%We do the Choleski decomposition because it is needed to compute random numbers from a multivariate gaussian distribution
    end
    for m = 1:M
        if (resample)
            start = prev_theta(:,m);             %if (resample), choose {X_{n-1}^i} accoding to weights {W_n^i}
            Eprev = LL(repl_ind(m));             %previous log likelihood(?)
            Pprev = max(sum(log(lognpdf(start,meanlogpar',sigmalogpar'))),-1e5); %post loglikelihood(?)
        else
            start = prev_theta(:,m);             %else choose the previous population {X_{n-1}^i} serially
            Eprev = LL(m);
            Pprev = max(sum(log(lognpdf(start,meanlogpar',sigmalogpar'))),-1e5);
        end
        x = zeros(parnum,burn);
        x(:,1) = start;
        acc = 1;
        while (acc<=burn)                       %MCMC
            prop = GMrnd_par(GMMtheta.mu,Chol_Sigma,GMMtheta.PComponents,1)';   %proposed point by the proposal distribution (gaussian mixture)
            propparams = exp(prop);
            [~, P_data] = ode45(@(t,x) toyModel(t,x,propparams(end-2),propparams(end-1),propparams(end),pathLengths,connectivity,propparams(1:end-3),[0,0],[0,0],[0,timeVector(end)]),timeVector,x0);

            Eprop = 0.5*sum((Y_data-P_data(:,observableNode)).^2)/sigmanoise;
            Pprop = max(sum(log(lognpdf(propparams,meanlogpar',sigmalogpar'))),-1e5);
            
            if rand < min(1,exp(beta*(Eprev-Eprop)-Pprev+Pprop)*GMpdf(GMMtheta,log(x(:,acc))')/GMpdf(GMMtheta,log(propparams)')*prod(propparams)/prod(x(:,acc)))
                x(:,acc+1) = propparams;
                Eprev = Eprop;
                Pprev = Pprop;
                accepts = accepts + 1;
            else
                x(:,acc+1) = x(:,acc);
            end
            acc = acc + 1;
        end
        next_theta(:,m) = x(:,end)';
        loglikelihood(b,m) = Eprev;
        incr_w(m) = exp(Eprev*(beta_schedule(b-1)-beta_schedule(b)));     %Why do we need this variable at all?
    end
    total_acc(b) = accepts;
    w(b+1,:) = incr_w;
    Theta{b} = next_theta;
    W(b+1,:) = W(b,:).*w(b+1,:)/sum(W(b,:).*w(b+1,:));             %Why do we multiply them by the previous weights?? calculatete weights {W_{n+1}^i}
    if (M/(1+var(W(b+1,:)/mean(W(b+1,:))))<0.5*M)                  %if they are degenerate, resample
        resample = 1;
        r_times = [r_times b+1];
        prob = W(b+1,:);
        W(b+1,:) = 1/M;
    else
        resample = 0;
    end
    minLL = min(loglikelihood(b,:));
end

%% Calculate evidence, gradient and Hessian
blength = length(beta_schedule);
M = size(Theta{end},2);

r71 = find(r_times==blength+1);
if ~isempty(r71)
    r_times = r_times(1:r71-1);
end
R = length(r_times)+1;
r_times = [1 r_times blength];
rterms = ones(R,M);
sterms = ones(1,R);
for k = 1:R
    for m = r_times(k)+1:r_times(k+1)
        rterms(k,:) = rterms(k,:).*w(m,:);
    end
    sterms(k) = mean(rterms(k,:));
end
marg_likelihood = prod(sterms);


ev = -marg_likelihood;
gF = ev*(mean(log(Theta{end}(variable_parameters,:)'))-mean(log(initial_points(variable_parameters,:)')));
Cprior = cov(log(initial_points(variable_parameters,:)'));
Cpost = cov(log(Theta{end}(variable_parameters,:)'));
H = ev*(Cpost-Cprior);
%%
% for i=1:30
%     for k = 1:length(meanlogpar)
%         subplot(length(meanlogpar),1,k)
%         hist(log(Theta{i}(k,:)),20);xlim([-6 6]);title(b)
%     end
%     pause(0.5)
% end


