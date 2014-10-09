clf

hyperparameters=-8:0.1:2;

likelihood=zeros(length(hyperparameters));%Likelihood computed just looking at a certain parameter value, and not using a prior distribution.
mu1=-2;
mu2=0;
sigma1=2;
sigma2=.1;
for i=1:length(hyperparameters)
    for j=1:length(hyperparameters)
        likelihood(i,j) = lognpdf(exp(hyperparameters(i)),mu1,sigma1)*lognpdf(exp(hyperparameters(j)),mu2,sigma2);
    end
end

imagesc(hyperparameters,hyperparameters,likelihood)
colorbar
set(gca,'YDir','normal')
xlabel('Strength2, logarithmic')
ylabel('Strength1, logarithmic')