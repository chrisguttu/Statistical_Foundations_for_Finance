function [pdf, cdf, themean, thevar, quant, theES]= ...
         GAt(xvec,d,v,theta, gammaforES)
% [pdf, cdf, themean, thevar, quant, theES] =GAt(xvec,d,v,theta, gammaforES);
% Here, quant is the value at risk (the gamma-quantile).
% The Generalized Asymmetric t, GAt, see page 273 of Intermediate
% Probability. 
% Observe that we need a scaled version to reduce to the Student's t.
% We need, for a given n>0, d=2, v=n/2 and theta=1, and
%  then if Z ~ GAt(d,v,theta), let X=sqrt(2)*Z.
% This program does NOT do the scaling!
% Note though: Calling:
%   n=4; dd=2; vv=n/2; 
%   [garb1, garb2, garb3, garb4, theES]=GAt(0,dd,vv,1,0.01); theES*sqrt(2)
% and then (for the regular Student's t): 
%   df=4; q=0.01; cc=tinv(q,df); ESt = -tpdf(cc,df) * (df+cc^2)/(df-1) / q
% yields the same value.

% Use this to check the Expected Shortfall calculation.
% d=2.0; v=1.5; theta=0.9; cutoff = 0.01;
% rr = GAtsim(50000, d,v,theta); 
% 
% [pdf, cdf, themean, thevar, theES] =GAt(0,d,v,theta, cutoff);
% q = GAtquantile(cutoff, d,v,theta);
% use=rr(rr<q); empirical01ES = mean(use), true01ES = theES
% 
% d=2.0; v=1.5; theta=0.9; cutoff = 0.05;
% [pdf, cdf, themean, thevar, theES] =GAt(0,d,v,theta, cutoff);
% q = GAtquantile(cutoff, d,v,theta);
% use=rr(rr<q); empirical05ES = mean(use), true05ES = theES

ll=length(xvec); xvec=reshape(xvec,ll,1); pdf = zeros(ll,1);
konst = 1 / ( beta(1/d,v) * v^(1/d) * (theta+1/theta) / d );
k = find(xvec<0);
if any(k), y=xvec(k); pdf(k) = ( 1 + (-y*theta).^d / v ).^(-(v+1/d)); end
k = find(xvec>=0);
if any(k), y=xvec(k); pdf(k) = ( 1 + (y/theta).^d / v ).^(-(v+1/d)); end
pdf = konst * pdf;
if nargout>1
  cdf=zeros(length(xvec),1);
  k = find(xvec<0);
  if any(k)
    y=xvec(k); L = v./(v+(-y*theta).^d);
    cdf(k) = betainc(L,v,1/d)/(1+theta^2);
  end
  k = find(xvec==0);
  if any(k), y=xvec(k); cdf(k) = 1/(1+theta^2); end
  k = find(xvec>0);
  if any(k)
    y=xvec(k); top=(y/theta).^d; U=top./(v+top);
    cdf(k) = 1/(1+theta^2) + betainc(U,1/d,v)/(1+theta^(-2));
  end
end


if nargout>2
   themean = NaN; thevar=NaN; theES=NaN;
   if v*d>1
     themean=GAtmom(1,d,v,theta); 
     quant = GAtquantile(gammaforES, d,v,theta);
     theES = GAttail(1,d,v,theta, quant); 
   end
   if v*d>2, thevar=GAtmom(2,d,v,theta) - themean^2; end
end

function m = GAtmom(r,d,v,theta)
t1 = ( (-1)^r * theta^(-r-1) + theta^(r+1) ) / (theta + 1/theta);
t2 = beta((r+1)/d , v-r/d) / beta(1/d,v) * v^(r/d);
m = t1 * t2;

function m = GAttail(r,d,v,theta, c)
t1 = (-1)^r * v^(r/d) * (1+theta^2) / (theta^r + theta^(r+2));
L = v / ( v+(-c*theta)^d );
t2 = beta(v-r/d, (r+1)/d) * betainc(L, v-r/d, ...
     (r+1)/d) / betainc(L,v,1/d) / beta(v,1/d);
m = t1*t2;
