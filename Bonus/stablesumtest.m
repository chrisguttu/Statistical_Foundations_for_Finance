function [ahat, tau, b90, b95, b99]=stablesumtest(z,sslen,doplot,Nperm)
% [ahat, tau, b90, b95, b99]=stablesumtest(z,sslen,doplot,Nperm)
% scale invariant test for null of location-zero iid symmetric stable data.
%
% IN: 
% z is a vector of (presumably) iid data of length between 500 and 10,000
%   SUCH THAT ITS LOCATION IS ZERO
%   (This can be done by using z=y-mu or z=(y-mu)/sigma, where mu and sigma are
%    estimated via McCulloch or MLE or other estimators.)
% sslen (optional) is the number of summed intervals to compute. Default given below.
% doplot (optional) is 1 to show a plot of the estimates.
% Nperm is number of permutations of the data set to use to compute the
%    test statistic, and then the average of them is used for conducting
%    the 3 cutoff tests. Default is 0.
%
% OUT: ahat is the Hint estimate of the stable tail index
%      tau is the value of the test statistic
%      b90 = 1 if we can reject the stable assumption at the 90% level, 0 otherwise
%      b95 = 1 "     "      "      "      "      "     "     95% level, 0 otherwise
%      b99 = 1 "     "      "      "      "      "     "     99% level, 0 otherwise
%
% Marc Paolella, January 1999

T=length(z);
ahat = Hint(z);
if T<500, warning('Inference on tau will be questionable for T much less than 1,500'), end
if T>10000, warning('Inference on tau will be questionable for T much greater than 10,000'), end

if nargin<2, sslen=-1; end
if nargin<3, doplot=0; end
if nargin<4, Nperm=0; end

if ahat >= 2 % then reject straight away
  ahat=2; tau=0; b90=1; b95=1; b99=1;
  return
end

if sslen<=0
  sslen=min(10,round(T/200));
  % override for 500<=T<=700
  if (T<700), sslen=5; end 
else
  if sslen<3, error ('sslen too small'), end  
  sslen=min(sslen,round(T/50));
end

[~,avec,stda,tau]=sumstab(z,sslen);
if Nperm>0
  ttau=zeros(Nperm,1);  
  for i=1:Nperm 
    ii=randperm(T); zperm=z(ii); [~, taup]=stablesumtest(zperm, sslen); ttau(i)=taup;  
  end
  tau=mean(ttau);
end
cutoff90=cutval(0.90,T,ahat); cutoff95=cutval(0.95,T,ahat); cutoff99=cutval(0.99,T,ahat);
b90=(tau>cutoff90); b95=(tau>cutoff95); b99=(tau>cutoff99);

if doplot==1
  xx=1:sslen;
  plot(xx,avec,'r-o',xx,avec+1.96*stda,'g--',xx,avec-1.96*stda,'g--','linewidth',2)
  grid,  ax=axis; axis([1 ax(2) ax(3) min(ax(4),2)])
  title (['\alpha = ',num2str(ahat),',  \tau = ',num2str(tau)]), set(gca,'Fontsize',18)
end  


function [ab,SE]=Hint(y)
% [ab,SE]=Hint(y)

nobs=length(y); nstar=nobs/1000;

step=max(round(nobs/100),1);
low=round(0.2*nobs); high=round(0.8*nobs);
kr=low:step:high; kl=length(kr); 
reg=(kr/1000)'; xmat=[ones(kl,1),reg];

h=hill(y,kr); beta=pinv(xmat'*xmat)*xmat'*h;
b=beta(1); m=beta(2);

q = [-0.8110,-0.3079,2.0278]; set=[1;b;sqrt(b)]; ab=q*set;

if nargout > 1
  SE= 0.03222739 -0.00204714 * nstar + 0.02273217 * (nstar^(-1)) ...
    -0.00083515 * (nstar^(-2));
end




function [len,a,stda,tstat]=sumstab(y,sumlen)
% [len,a,stda,tstat]=sumstab(y,estim,sumlen)
% default sumlen is such that the last value has 100 observations.

nobs=length(y);
if nargin<2
  sumlen=round(nobs/100);
end
if length(sumlen)==1, sumlen=1:sumlen; end

len=0; a=0; b=0; stda=0; stdb=0;
[r,c]=size(y); if c>1, y=y'; end
for i=1:length(sumlen)
  n=sumlen(i);
  if 1==2 % the elegant, but memory hungry way
    set=y(1:n*floor(nobs/n));
    nice=kron(eye(floor(nobs/n)),ones(1,n));
    rn=nice * set;
    len(i)=length(rn);
  else
    len(i) = fix(nobs/n);
    rn=zeros(len(i),1);
    for j=1:len(i)
      rn(j)=sum(y((j-1)*n+1 : n*j));
    end
  end
  [a(i),stda(i)]=Hint(rn);
end

if nargout>2
  Oi=diag(stda.^(-2)); X=[ones(length(a),1), sumlen'];
  reg=inv(X'*Oi*X)*X'*Oi*a';
  C=diag(1./stda); 
  ehat=C*(a' - X*reg); 
  shat=(1/(length(a)-2)) * sum(ehat.^2);
  varcov=shat*inv(X'*Oi*X);

  se=sqrt(shat); 
  tstat=reg(2)/sqrt(varcov(2,2));
end




function C=cutval(gamma,T,alpha)
% C=cutval(gamma,T,alpha), alpha can be a vector.

T2=sqrt(T);
if (gamma==0.90)

  c0 =  6.89125873970754  -0.00022591029043*T  +0.00373170040702*T2;
  c1 = -3.40461054740374  +0.00008045493031*T  +0.04676122520617*T2;
  c2 =  0.49863480029558  +0.00005479817898*T  -0.02744839317119*T2;
  
elseif (gamma==0.95)

  c0 =  8.37678998539416  -0.00010067103790*T  -0.00161588975535*T2;
  c1 = -1.40754883873540  +0.00038618206953*T  +0.00136231416456*T2;
  c2 = -1.12358924203789  -0.00023341297751*T  +0.01119915597888*T2;

elseif (gamma==0.99)

  c0 = 14.27777312550145  -0.00002561601571*T  +0.00778459157178*T2;
  c1 = -4.08083390855784  +0.00002768864957*T  +0.02061453596457*T2;
  c2 = -0.58387014348690  -0.00000548302271*T  -0.01109110722591*T2;
                                             
else
  error ('bad gamma level')
end

C = c0 + c1*alpha + c2*alpha.^2;


function [hill,hstd]=hill(x,krange)
% hill=hill(x , krange=see below)
% The Hill tail estimate.

% n=length(x);  

y=sort(abs(x(x~=0))); n=length(y);
lny=log(y); % just for the Hill estimator

for loop=1:length(krange) 
  k=krange(loop);
  a=1/(  mean(lny(  (n+1-k):n  )) - lny(n-k)  );
  hillmat(loop)= a;
  vmat(loop)=a*a/k;
end
hill=hillmat';
hstd=sqrt(vmat');
