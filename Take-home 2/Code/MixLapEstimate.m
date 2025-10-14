function [param,stderr,iters,loglik,Varcov] = MixLapEstimate(x, initvec, lmda)
% param: (mu11, mu21, mu12, mu22,
%         S1_11, S1_12, S1_22, S2_11, S2_12, S2_22, 
%         b1, b2)
[nobs, d]=size(x); 
if d~=2 
    error('Not done yet, use EM instead.')
end

if d==2
% mu11, mu21, mu21, mu22,
% S1_11, S1_12, S1_22, S2_11, S2_12, S2_22, 
% b1, b2
bound.lo = [-1 -1 -1 -1 ...
            0.01 -90 0.01 0.01 -90 0.01 ...
            0.01 0.01];
bound.hi = [1 1 1 1 ...
            90 90 90 90 90 90 ...
            20 20];
bound.which = [0 0 0 0 ...
               1 1 1 1 1 1 ...
               1 1];

    if nargin == 1
        initvec =[-0.8 -0.2 -0.8 -0.2 20 2 10 20 2 10 8 4];
        lmda = 0.5;
    end
end

maxiter=500; tol=1e-8;
MaxFunEvals = length(initvec)*maxiter;
opts = optimset('Display', 'iter', 'Maxiter', maxiter, 'TolFun', tol, ...
    'TolX', tol, 'MaxFunEvals', MaxFunEvals, 'LargeScale', 'Off');

[pout, fval, ~, theoutput, ~, hess] = ...
fminunc(@(param) MixLaploglik(param, x, lmda, bound), ...
einschrk(initvec, bound), opts);
V = inv(hess)/nobs; % Don't negate because we work with the negative of the loglik
[param, V] = einschrk(pout, bound, V); % transform and apply delta method to get V
param = param'; 
Varcov = V; 
stderr = sqrt(diag(V)); % Approximate standard errors
loglik = -fval*nobs; 
iters = theoutput.iterations;

function ll = MixLaploglik(param, x, lmda, bound)
if nargin<4, bound=0; end
if isstruct(bound), param=einschrk(real(param),bound,999); end
[nobs, d] = size(x); 
Sig1 = zeros(d,d);
Sig2 = zeros(d,d);
mu1 = param(1:2); % Assume d=2
mu2 = param(3:4);

Sig1(1,1)=param(5); Sig1(1,2)=param(6); 
Sig1(2,2)=param(7); Sig1(2,1)=Sig1(1,2);

Sig2(1,1)=param(8); Sig2(1,2)=param(9); 
Sig2(2,2)=param(10); Sig2(2,1)=Sig2(1,2);

b1 = param(11);
b2 = param(12);

if min([eig(Sig1), eig(Sig2)], [], "all") < 1e-10, ll = 1e5;
else
pdf = zeros(nobs,1);
for i = 1:nobs
    pdf(i) = lmda*mvlpdf(x(i,:), mu1, Sig1, b1) + ... 
    (1-lmda)*mvlpdf(x(i,:), mu2, Sig2, b2); 
end
llvec=log(pdf); ll=-mean(llvec); if isinf(ll), ll=1e5; end
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

function y = mvlpdf(x, mu, Sigma, b)
% Joint PDF of multivariate Laplace based on III. (14.31) 
d = length(x);
x = reshape(x,d,1); 
mu = reshape(mu,d,1); 
m = (x-mu)' * (Sigma \(x-mu));
y = exp(log(2) + (b/2 - d/4) * log(m/2) + log(besselk(b-d/2, sqrt(2*m))) - ...
    log(det(Sigma))/2 - log(2*pi)*d/2 - gammaln(b));