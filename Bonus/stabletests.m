function [Ap, tp, Cp, At, tt]=stabletests(x,Asym,Nperm)
% [Ap, tp, Cp, At, tt]=stabletests(x,Asym=0,Nperm=1e3)
%
% INPUT:
%   x is data vector (of length T)
%   Asym=0 (default) to assume symmetry.
%     otherwise, performs invcdf transform
%     and requires path to Nolan stable toolbox.
%   Nperm is number of permutations to be used in the tau test.
%
% OUTPUT
% Ap:   Alhadi p-value 
% tp:   tau_{Nperm} p-value 
% Cp:  A+tau_{Nperm} p-value
% At: Alhadi test statistic
% tt: tau_{Nperm} test statistic
%
%
%   Alhadi      test is availble for 250<=T<=10000
%   tau_{Nperm} test is availble for 500<=T<=10000
%   Combined    test is availble for T=500 and T=1000.
%     (For other sample sizes for Combined, need to simulate)
%
% Marc S. Paolella, August 2015.
%
% Example. For size check.
%   cut=0.05; alpha=1.7; beta=0; T=1000; sim=1e3; 
%   AAp=zeros(sim,1); ttp=AAp; CCp=AAp;
%   for i=1:sim , if mod(i,10)==0, i, end
%     x=stabgen(T, alpha,beta,1,0,i);
%     [Ap, tp, Cp]=stabletests(x,0,20);
%     AAp(i)=Ap; ttp(i)=tp; CCp(i)=Cp;
%   end
%   size05=[mean(AAp<cut), mean(ttp<cut), mean(CCp<cut)] 
%   
% Example. For power against Student's t
%   cut=0.05; df=4; T=1000; sim=1e3; 
%   AAp=zeros(sim,1); ttp=AAp; CCp=AAp;
%   for i=1:sim , if mod(i,10)==0, i, end
%     x=trnd(df,[T 1]);
%     [Ap, tp, Cp]=stabletests(x,0,20);
%     AAp(i)=Ap; ttp(i)=tp; CCp(i)=Cp;
%   end
%   size05=[mean(AAp<cut), mean(ttp<cut), mean(CCp<cut)] 
  
if nargin<2, Asym=0; end
Ap=NaN; tp=NaN; Cp=NaN; At=NaN; tt=NaN;  %#ok<NASGU>

T=length(x); 
if (T==500),  load COMBTESTAlhaditauT500, end
if (T==1000), load COMBTESTAlhaditauT1000, end

if Asym~=0 
  theta = stablefit(x,1,1); % 2nd param is method: 1 is MLE, 3 is cf
  alpha_hat=theta(1); beta_hat=theta(2); sigma_hat=theta(3); mu_hat=theta(4);
  x2 = (x - mu_hat) / sigma_hat; 
  F = stableqkcdf(x2,[alpha_hat,beta_hat,1,0],1);
  usex = stableqkinv(F,[alpha_hat,0,1,0],1);
else
  usex=x;
end
[~, tt, tp]=stablesumtestwithpvalue(usex,Nperm);
[At, Ap]=stableALHADItestwithpvalue(usex);
if (((T==500) || (T==1000)) && ((Ap>0) && (tp>0)))
  pC=log(Ap*tp); pCvec=zeros(length(alphavec),1); 
  for j=1:length(alphavec), pCvec(j)=mean(Cvec(j,:)<=pC); end %#ok<NODEF>
  aHint=Hint(usex); Cp=interp1(alphavec,pCvec,aHint);
end
