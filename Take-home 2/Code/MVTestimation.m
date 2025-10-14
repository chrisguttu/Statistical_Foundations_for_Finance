function [param,stderr,iters,loglik,Varcov] = MVTestimation(x, initvec)
% param: (k, mu1, mu2, Sigma_11, Sigma_12, Sigma_22)
[nobs, d]=size(x); if d~=2, error('not done yet, use EM'), end
if d==2
%%%%%%%% k mu1 mu2 s11 s12 s22
bound.lo= [ 0.2 -1 -1 0.01 -90 0.01];
bound.hi= [ 20 1 1 90 90 90];
bound.which=[ 1 0 0 1 1 1];
if nargin < 2
initvec =[2 -0.8 -0.2 20 2 10];
end
end
maxiter=300; tol=1e-7; MaxFunEvals=length(initvec)*maxiter;
opts=optimset('Display','iter','Maxiter',maxiter,'TolFun',tol,'TolX',tol,...
'MaxFunEvals',MaxFunEvals,'LargeScale','Off');
[pout,fval,~,theoutput,~,hess]= ...
fminunc(@(param) MVTloglik(param,x,bound),einschrk(initvec,bound),opts);
V=inv(hess)/nobs; % Don't negate because we work with the negative of the loglik
[param,V]=einschrk(pout,bound,V); % transform and apply delta method to get V
param=param'; Varcov=V; stderr=sqrt(diag(V)); % Approximate standard errors
loglik=-fval*nobs; iters=theoutput.iterations;

function ll=MVTloglik(param,x,bound)
if nargin<3, bound=0; end
if isstruct(bound), param=einschrk(real(param),bound,999); end
[nobs, d]=size(x); Sig=zeros(d,d); k=param(1); mu=param(2:3); % Assume d=2
Sig(1,1)=param(4); Sig(2,2)=param(6); Sig(1,2)=param(5); Sig(2,1)=Sig(1,2);
if min(eig(Sig))<1e-10, ll=1e5;
else
pdf=zeros(nobs,1);
for i=1:nobs, pdf(i) = mvtpdfmine(x(i,:),k,mu,Sig); end
llvec=log(pdf); ll=-mean(llvec); if isinf(ll), ll=1e5; end
end

function [pout, Vout] = einschrk(pin, bound, Vin)
lo = bound.lo; hi = bound.hi ; welche = bound.which ;
if nargin < 3
trans=sqrt((hi-pin) ./ (pin-lo)); pout=(1-welche) .* pin + welche .* trans ;
Vout =[];
else
trans=(hi+lo .* pin.^2) ./ (1 + pin.^2); pout=(1-welche) .* pin + welche .* trans;
% now adjust the standard errors
trans=2*pin .* (lo-hi) ./ (1+pin.^2).^2;
d=(1-welche) + welche .* trans; % either unity or delta method.
J=diag(d); Vout = J * Vin * J;
end

function y = mvtpdfmine(x,df,mu,Sigma)
% x is a d X 1 vector. Unlike Matlab's version, cannot pass a matrix.
% Matlab's routine accepts a correlation (not dispersion) matrix.
% So, just need to do the usual scale transform. For example:
% x=[0.2 0.3]'; C = [1 .4; .4 1]; df = 2;
% scalevec=[1 2]'; xx=x./scalevec; mvtpdf(xx,C,df)/prod(scalevec)
% Same as:
% Sigma = diag(scalevec) * C * diag(scalevec); mvtpdfmine(x,df,[],Sigma)
d=length(x);
if nargin<3, mu = []; end, if isempty(mu), mu = zeros(d,1); end
if nargin<4, Sigma = eye(d); end
x = reshape(x,d,1); mu = reshape(mu,d,1); term = (x-mu)' * inv(Sigma) * (x-mu);
logN=-((df+d)/2)*log(1+term/df); logD=0.5*log(det(Sigma))+(d/2)*log(df*pi);
y = exp(gammaln((df+d)/2) - gammaln(df/2) + logN - logD);