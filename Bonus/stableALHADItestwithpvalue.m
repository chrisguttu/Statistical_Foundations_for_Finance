function [alhaditeststat, pval]=stableALHADItestwithpvalue(x,useahat)
% [alhaditeststat, pval]=stableALHADItestwithpvalue(x,useahat)
%
% For a given vector data set x (with 250<=N<=20000) 
%   and alpha-hat in [1.05, 1.99]
% Get the test statistic of the ALHADI stable test. 
% And then for all N,alpha combinations, compute the p-value.
% Then apply double interpolation to get the delivered p-value.
%
% do not pass fixahat: It is just for assessing the size if somehow
%   the true alpha is known. See the paper.

T=length(x); 
if (T<250) || (T>20000), error('Sample size out of range'), end 

[~,~,sigma_hat,mu_hat]=stablecull(x); x = (x - mu_hat) / sigma_hat; 
ahat=Hint(x); alhaditeststat = ahat-stablecull(x);

% fix the value of alpha, for an exercise in the paper for assessing the size
if nargin<2, useahat=ahat; end

if ahat<1.05, pval=-1;  return, end
if ahat>1.99, pval=-1;  return, end 

if 1==2
  load ('MixNparammatrixA', 'MixNparammatrixA', 'Tvec', 'Tlen', 'alphavec', 'alphalen')
  PMAT=MixNparammatrixA;
else
  load ('MixNparammatrixANEW', 'MixNparammatrixANEW', 'Tvec', 'alphavecNEW')
  PMAT=MixNparammatrixANEW; alphavec=alphavecNEW;
end
alphalen=length(alphavec);
if ismember(T,Tvec) % Is the sample size precisely one in the table?
  pos=find(T==Tvec); pvalMat=zeros(1,alphalen);
  for aa=1:alphalen
    param=reshape(PMAT(pos,aa,:), 5,1); %#ok<FNDSB>
    mu1=param(1); mu2=param(2); sig1=param(3); sig2=param(4); lam=param(5);
    z1=(alhaditeststat-mu1)/sig1; cdf1=normcdf(z1);
    z2=(alhaditeststat-mu2)/sig2; cdf2=normcdf(z2);
    pvalMat(aa) = 1 - (lam*cdf1+(1-lam)*cdf2);
  end
  pval=interp1(alphavec,pvalMat,useahat);
else
  Tlen=length(Tvec);
  pvalMat=zeros(Tlen,alphalen);
  for TT=1:Tlen
    for aa=1:alphalen
      param=reshape(PMAT(TT,aa,:), 5,1);
      % just to look. [Tvec(TT) alphavec(aa) param']
      mu1=param(1); mu2=param(2); sig1=param(3); sig2=param(4); lam=param(5);
      z1=(alhaditeststat-mu1)/sig1; cdf1=normcdf(z1);
      z2=(alhaditeststat-mu2)/sig2; cdf2=normcdf(z2);
      pvalMat(TT,aa) = 1 - (lam*cdf1+(1-lam)*cdf2);
    end
  end
  pval=interp2(Tvec,alphavec,pvalMat',T,useahat);
end
