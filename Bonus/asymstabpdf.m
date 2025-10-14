function f = asymstabpdf(xvec, a, b, plotintegrand)
% pdf of the asymmetric stable. See also asymstab.m
% Set plotintegrand to 1, and xvec a scalar, to plot the integrand.

if nargin<4, plotintegrand=0; end
if a==1, error('not ready yet'), end

usevectorized =1; % vectorized is faster, but not massively
xl = length(xvec); f = zeros(xl, 1);

if usevectorized
  k=find(xvec==0);  
  if ~isempty(k), pdf0 = stab(0, a, b); f(k)=pdf0; end
  k=find(xvec>0);
  if ~isempty(k), f(k)=stabvec(xvec(k),a,b); end
  k=find(xvec<0);
  if ~isempty(k), f(k)=stabvec(-xvec(k),a,-b); end
else
  for loop = 1:xl
    x = xvec(loop);
    if (x < 0), f(loop) = stab(-x, a, -b); else f(loop) = stab(x, a, b); end
  end
end

if plotintegrand && xl==1 % show a plot of the integrand
  if x<0, x=-x; b=-b; end
  t0 = (1 / a) * atan(b * tan((pi * a) / 2));
  t=-t0:0.001:(pi/2); I=integrand2(t, x, a, t0); plot(t,I,'r-')
end

function pdf = stabvec(xvec, a, b)
t0 = (1 / a) * atan(b * tan((pi * a) / 2));
if (xvec == 0) % SFB 649 Discussion Paper 2005-008
  pdf = gamma(1 + (1 / a)) * cos(t0) / pi / (1 + (b * tan((pi * a) / 2))^2)^(1/ (2 * a));
else % from Stoyan et.al.
  xvecabs=abs(xvec);
  tol = 1e-9; display = 0;
  integ = integral(@(t) integrand(t, xvecabs, a, t0), -t0, pi/2,'ArrayValued',true);
  pdf = a * xvec.^(1 / (a-1)) / pi / abs(a-1) .* integ;
end

function pdf = stab(x, a, b)
t0 = (1 / a) * atan(b * tan((pi * a) / 2));
if (x == 0) % SFB 649 Discussion Paper 2005-008
  pdf = gamma(1 + (1 / a)) * cos(t0) / pi / (1 + (b * tan((pi * a) / 2))^2)^(1/ (2 * a));
else % from Stoyan et.al.
  tol = 1e-9; display = 0;
  integ = quadv(@integrand, -t0, pi/2, tol, display, abs(x), a, t0);
  pdf = a * x.^(1 / (a-1)) / pi / abs(a-1) .* integ;
end

function I = integrand(t, x, a, t0)
ct = cos(t); s = t0 + t;
v = (cos(a * t0))^(1 / (a-1)) * (ct / sin(a*s))^(a / (a-1)) * (cos(a*s - t) / ct);
term = -x.^(a / (a-1)); I = v .* exp(term * v);

function I = integrand2(t, x, a, t0)
ct = cos(t); s = t0 + t;
v = (cos(a * t0))^(1 / (a-1)) * (ct / sin(a*s)).^(a / (a-1)) .* (cos(a*s - t) ./ ct);
term = -x.^(a / (a-1)); I = v .* exp(term * v);
