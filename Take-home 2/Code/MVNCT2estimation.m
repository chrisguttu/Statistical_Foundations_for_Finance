function [param,stderr,iters,loglik,Varcov] = MVNCT2estimation(x, initvec)
[d, T]=size(x'); if d~=2, error('not done yet, use 2-step'), end
%%%%%%%% k mu1 mu2 scale1 scale2 R12 gam1 gam2
bound.lo = [ 1 -10 -10 0.01 0.01 -1 -4 -4 ];
bound.hi = [ 20 10 10 90 90 1 4 4 ];
bound.which = [ 1 1 1 1 1 1 1 1 ];
if nargin < 2
    initvec = [ 3 0 0 2 2 0.5 0 0 ];
end
maxiter=500; tol=1e-10; MaxFunEvals=length(initvec)*maxiter;
opts=optimset('Display','iter','Maxiter',maxiter,'TolFun',tol,'TolX',tol,...
'MaxFunEvals',MaxFunEvals,'LargeScale','Off');
[pout,fval,~,theoutput,~,hess]= ...
fminunc(@(param) MVNCTloglik(param,x',bound),einschrk(initvec,bound),opts);
V=inv(hess)/T; [param,V]=einschrk(pout,bound,V); param=param';
Varcov=V; stderr=sqrt(diag(V)); loglik=-fval*T; iters=theoutput.iterations;

function ll=MVNCTloglik(param,x,bound)
if nargin<3, bound=0; end
if isstruct(bound), param=einschrk(real(param),bound,999); end
k=param(1); mu=param(2:3); scale=param(4:5); gam=param(7:8);
R12=param(6); R=[1 R12; R12 1]; 

if min(eig(R))<1e-10, ll=1e5;
else
xx=x; for i=1:2, xx(i,:)=(x(i,:)-mu(i))/scale(i); end
mu_tr = [0, 0];
llvec = mvnctpdfln(xx, mu_tr, gam, k, R) - log(prod(scale));
ll=-mean(llvec); if isinf(ll), ll=1e5; end
end

function [pout, Vout] = einschrk(pin, bound, Vin)
lo = bound.lo; hi = bound.hi ; welche = bound.which;
if nargin < 3
trans=sqrt((hi-pin) ./ (pin-lo)); pout=(1-welche) .* pin + welche .* trans ;
Vout =[];
else
trans=(hi+lo .* pin.^2) ./ (1 + pin.^2); pout=(1-welche) .* pin + welche .* trans;
% now adjust the standard errors
trans=2*pin .* (lo-hi) ./ (1+pin.^2).^2;
d=(1-welche) + welche .* trans; % either unity or delta method.
J=diag(d); Vout = J*Vin*J;
end

function pdfln = mvnctpdfln(x, mu, gam, v, Sigma)
% x d X T matrix of evaluation points
% mu, gam d-length location and noncentrality vector
% v is df; Sigma is the dispersion matrix.
[d,t] = size(x); C=Sigma; [R, err] = cholcov(C, 0);
assert(err == 0, 'C is not (semi) positive definite');
mu=reshape(mu,length(mu),1); gam=reshape(gam,length(gam),1);
vn2 = (v + d) / 2; xm = x - repmat(mu,1,t); rho = sum((R'\xm).^2,1);
pdfln = gammaln(vn2) - d/2*log(pi*v) - gammaln(v/2) - ...
sum(slog(diag(R))) - vn2*log1p(rho/v);
if (all(gam == 0)), return; end
idx = (pdfln >= -37); maxiter=1e4; k=0;
if (any(idx))
gcg = sum((R'\gam).^2); pdfln = pdfln - 0.5*gcg; xcg = xm' * (C \ gam);
term = 0.5*log(2) + log(xcg) - 0.5*slog(v+rho');
term(term == -inf) = log(realmin); term(term == +inf) = log(realmax);
logterms = gammaln((v+d+k)/2) - gammaln(k+1) - gammaln(vn2) + k*term;
ff = real(exp(logterms)); logsumk = log(ff);
while (k < maxiter)
k=k+1;
logterms = gammaln((v+d+k)/2) - gammaln(k+1) - gammaln(vn2) + k*term(idx);
ff = real(exp(logterms-logsumk(idx))); logsumk(idx)=logsumk(idx)+log1p(ff);
idx(idx) = (abs(ff) > 1e-4); if (all(idx == false)), break, end
end
pdfln = real(pdfln+logsumk');
end

function y = slog(x) % Truncated log. No -Inf or +Inf.
y = log(max(realmin, min(realmax, x)));