 function   [log_lh,mahalaD]=wdensity(X, mu, Sigma, p, sharedCov, CovType)
 %WDENSITY Weighted conditional density and mahalanobis distance.
 %   LOG_LH = WDENSITY(...) returns log of component conditional density
 %   (weighted by the component probability) of X. LOG_LH is a N-by-K matrix
 %   LOG_LH, where K is the number of Gaussian components. LOG_LH(I,J) is
 %   log (Pr(point I|component J) * Prob( component J))
 %
 %   [LOG_LH, MAHALAD]=WDENSITY(...) returns the Mahalanobis distance in
 %   the N-by-K matrix MAHALAD. MAHALAD(I,J) is the Mahalanobis distance of
 %   point I from the mean of component J. 
 %   Copyright 2007 The MathWorks, Inc.
 %   $Revision: 1.1.6.3 $  $Date: 2008/10/08 18:23:20 $
 
     log_prior = log(p);
     [n,d]=size(X);
     k=size(mu,1);
     log_lh = zeros(n,k);
     mahalaD = zeros(n,k);
     logDetSigma = -Inf;
     for j = 1:k
         if sharedCov
             if j == 1
                 if CovType == 2 % full covariance
                     [L,f] = chol(Sigma);
                     diagL = diag(L);
                     if (f ~= 0)|| any(abs(diagL) < eps(max(abs(diagL)))*size(L,1))
                         error('stats:gmdistribution:wdensity:IllCondCov', ...
                             'Ill-conditioned covariance created.');
                     end
                     logDetSigma = 2*sum(log(diagL));
                 else %diagonal
                     L = sqrt(Sigma);
                     if  any(L < eps( max(L))*d)
                         error('stats:gmdistribution:wdensity:IllCondCov', ...
                             'Ill-conditioned covariance created.');
                     end
                     logDetSigma = sum( log(Sigma) );
                 end
             end
         else %different covariance
             if CovType == 2 %full covariacne
                 % compute the log determinant of covariance
                 [L,f] = chol(Sigma(:,:,j) );
                 diagL = diag(L);
%                  if (f ~= 0) || any(abs(diagL) < eps(max(abs(diagL)))*size(L,1))
%                      error('stats:gmdistribution:wdensity:IllCondCov', ...
%                          'Ill-conditioned covariance created');
%                  end
                 logDetSigma = 2*sum(log(diagL));
             else %diagonal covariance
                 L = sqrt(Sigma(:,:,j)); % a vector
                 if  any(L < eps(max(L))*d)
                     error('stats:gmdistribution:wdensity:IllCondCov', ...
                         'Ill-conditioned covariance created.');
                 end
                 logDetSigma = sum( log(Sigma(:,:,j)) );
 
             end
         end
 
         %Xcentered = bsxfun(@minus, X, mu(j,:));
         Xcentered = X-mu(j*ones(n,1),:);
         if CovType == 2
             xRinv = Xcentered /L ;
             %xRinv = solve_tril(L',Xcentered')';
         else
             xRinv = bsxfun(@times,Xcentered , (1./ L));
         end
         cols = size(xRinv,2);
         xRinv2 = xRinv.^2;
         mahalaD(:,j) = xRinv2*ones(cols,1);
         %mahalaD(:,j) = sum(xRinv.^2, 2);
 
         log_lh(:,j) = -0.5 * mahalaD(:,j) +...
             (-0.5 *logDetSigma + log_prior(j)) - d*log(2*pi)/2;
         %get the loglikelihood for each point with each component
         %log_lh is a N by K matrix, log_lh(i,j) is log \alpha_j(x_i|\theta_j)
     end
  