function pnew=condp(pin)

p = pin./max(pin(:));
p = p+eps;
pnew=p./repmat(sum(p,1),size(p,1),1);