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


function [hill,hstd]=hill(x,krange)
% hill=hill(x , krange=see below)
% The Hill tail estimate.

% n=length(x);  

%use the absolute value because we assume the distribution to be symmetric
y=sort(abs(x(x~=0))); n=length(y);
lny=log(y); % just for the Hill estimator
hillmat=zeros(1,length(krange)); vmat=hillmat;
for loop=1:length(krange) 
  k=krange(loop);
  a=1/(  mean(lny(  (n+1-k):n  )) - lny(n-k)  );
  hillmat(loop)= a; vmat(loop)=a*a/k;
end
hill=hillmat'; hstd=sqrt(vmat');

