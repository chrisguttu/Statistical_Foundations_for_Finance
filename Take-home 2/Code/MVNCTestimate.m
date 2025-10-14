function [param,stderr,iters,loglik,Varcov] = MVNCTestimate(x, initvec)
% param: (mu1, mu2, gam1, gam2, v, S11, S12, S22)
[nobs, d]=size(x); 
if d~=2 
    error('Not done yet, use EM instead.')
end

if d==2
% mu1, mu2, gam1, gam2, v, S11, S12, S22
bound.lo = [-10 -10 -4 -4 1 0.01 -90 0.01];
bound.hi = [10 10 4 4 20 90 90 90];
bound.which = [1 1 1 1 1 1 1 1];

if nargin < 2
    initvec =[0.1 0.1 0.1 0.8 3 10 2 10];
end
end

maxiter=500; tol=1e-10;
MaxFunEvals = length(initvec)*maxiter;
opts = optimset('Display', 'iter', 'Maxiter', maxiter, 'TolFun', tol, ...
    'TolX', tol, 'MaxFunEvals', MaxFunEvals, 'LargeScale', 'Off');

[pout, fval, ~, theoutput, ~, hess] = ...
fminunc(@(param) MVNCTloglik(param, x, bound), ...
einschrk(initvec, bound), opts);
V = inv(hess)/nobs; % Don't negate because we work with the negative of the loglik
[param, V] = einschrk(pout, bound, V); % transform and apply delta method to get V
param = param'; 
Varcov = V; 
stderr = sqrt(diag(V)); % Approximate standard errors
loglik = -fval*nobs; 
iters = theoutput.iterations;

function ll = MVNCTloglik(param, x, bound)
if nargin<3, bound=0; end
if isstruct(bound), param=einschrk(real(param),bound,999); end
d = length(x(1,:)); 
Sig = zeros(d,d);
mu = param(1:2); % Assume d=2
gam = param(3:4);
v = param(5);
Sig(1,1)=param(6); Sig(1,2)=param(7); 
Sig(2,2)=param(8); Sig(2,1)=Sig(1,2);

if min(eig(Sig)) < 1e-10, ll = 1e5;
else
llvec = mvnctpdfln(x', mu, gam, v, Sig); 
ll=-mean(llvec); 
if isinf(ll), ll=1e5; end
end

function [pout, Vout] = einschrk(pin, bound, Vin)
lo = bound.lo; hi = bound.hi; welche = bound.which;
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