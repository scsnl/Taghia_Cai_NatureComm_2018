function anew=logsumexp(a,b)
amax=max(a); 
A =size(a,1);
anew = amax + log(sum(exp(a-repmat(amax,A,1)).*b)+1e-310);