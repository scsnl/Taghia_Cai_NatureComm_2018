
n = size(Y,2);
p = size(Y,1);
s = size(Lm,2);
nSubjs = length(Ycell);
for ns = 1:nSubjs
col = 0; logQns = [];
for t = 1:s
  col = col + 1;
  kt = size(Lm{t},2);
  LmpsiiLm = Lm{t}'*diag(psii+eps)*Lm{t};
  temp = LmpsiiLm + reshape(reshape(Lcov{t},kt*kt,p)*psii,kt,kt);
  tempXm = Xm{t}(:,1+(ns-1)*n/nSubjs:(ns*n/nSubjs));
  logQns(:,col) = -.5*( +sum(Ycell{ns}.*(diag(psii+eps)*(Ycell{ns}-2*Lm{t}*tempXm)),1)' ...
      +reshape(temp,kt*kt,1)'*reshape(Xcov{t},kt*kt,1) ...
      +sum( tempXm.*(temp*tempXm) ,1)' ...
      +trace(Xcov{t}(2:end,2:end)) ...
      +sum( tempXm(2:end,:).*tempXm(2:end,:) ,1)' ...
      -2*sum(log(diag(chol(Xcov{t}(2:end,2:end)+eps)))) ...
      );
end %t
logOutProbs{ns} = logQns';
end

