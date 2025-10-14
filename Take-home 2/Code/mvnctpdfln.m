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