%res = digamma(X)
%
% Calculates the digamma function.
%
% Loops on rows of a matrix - could be optimized better
% 
% Matthew J. Beal GCNU 06/02/01

% This is a hack (Rasmussen).

function res=digamma(X);

coef=[-1/12 1/120 -1/252 1/240 -1/132 691/32760 -1/12];
krange = [1:7]';


[n, m]=size(X);
res=zeros(n,m);

for i=1:n,
  x=X(i,:);
  y = ceil(max(0,6-x));
  
  z = x + y;
  
  logz = log(z);
  
  res(i,:) = logz - .5./z + coef*exp(-2*krange*logz) ...
      - (y>=1)./x - (y>=2)./(x+1) - (y>=3)./(x+2) - (y>=4)./(x+3) - ...
      (y>=5)./(x+4) - (y>=6)./(x+5);
end;
  

