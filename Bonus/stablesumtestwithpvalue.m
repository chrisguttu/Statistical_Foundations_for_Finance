function [ahat, tau, pval]=stablesumtestwithpvalue(x,Nperm)
% [ahat, tau, pval]=stablesumtestwithpvalue(x,Nperm=20)
%
% For a given vector data set x (with 1000<=N<=10000) 
%   and alpha-hat in [1.1, 1.99]
% Get the test statistic of the summability stable test. 
% And then for all N,alpha combinations, compute the p-value.
% Then apply double interpolation to get the delivered p-value.

if nargin<2, Nperm=20; end
N=length(x); 
if (N<500) || (N>10000), error('Sample size out of range'), end 

[~,~,sigma_hat,mu_hat]=stablecull(x); x = (x - mu_hat) / sigma_hat; 

[ahat, tau]=stablesumtest(x,0,0,Nperm);
if ahat<1.1, ahat=1.1; end 
if ahat>1.99, ahat=1.99; end 
load MixGAtparammatrixtau
MixGAtparammatrixtau=MixGAtparammatrixtauNEW;
Nlen=length(Nvec); alphalen=length(alphavec);

pvalMat=zeros(Nlen,alphalen);
for NN=1:Nlen
  for aa=1:alphalen
    param=reshape(MixGAtparammatrixtau(NN,aa,:), 11,1);
    d1=param(1); v1=param(2); theta1=param(3); mu1=param(4); c1=param(5);
    d2=param(6); v2=param(7); theta2=param(8); mu2=param(9); c2=param(10);
    lam1=param(11);
    %theN=Nvec(NN), thealpha=alphavec(aa), d1, v1, theta1
    z1=(tau-mu1)/c1; [~, cdf1]= GAt(z1,d1,v1,theta1); 
    z2=(tau-mu2)/c2; [~, cdf2]= GAt(z2,d2,v2,theta2);
    pvalMat(NN,aa) = 1 - (lam1*cdf1+(1-lam1)*cdf2);
  end
end
pval=interp2(Nvec,alphavec,pvalMat',N,ahat);

