function logbeta1 = VBHMMbackward(v,phghm,logOutput)

T = size(v,2); 
H = size(logOutput,1); 
logbeta1 = zeros(H,T);
logbeta1(:,T) = zeros(H,1);
for t=T:-1:2
    phatvgh1 = zeros(H,1);
   for h = 1:H
        phatvgh1(h) = exp(logOutput(h,t))+1e-100; % eps
   end
    logbeta1(:,t-1)=logsumexp(repmat(logbeta1(:,t),1,H),repmat(phatvgh1,1,H).*phghm);
end
if any(isnan(logbeta1(:)))
        error('logbeta1= NaN');
end;

if any(isinf(logbeta1(:)))
        error('logbeta1= Inf');
end;
