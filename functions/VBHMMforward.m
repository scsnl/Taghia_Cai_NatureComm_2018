function [logalpha1,loglik1] = VBHMMforward(v,phghm,ph1,logOutput)

T = size(v,2); 
H = size(logOutput,1);  
logalpha1 = zeros(H,T);
for h = 1:H
    logalpha1(h,1) = logOutput(h,1) + log(ph1(h)) +1e-100; % eps
end
for t = 2:T
    phatvgh1 = zeros(H,1);
    for h = 1:H
        phatvgh1(h) = exp(logOutput(h,t)) + 1e-100; % eps 
    end
    logalpha1(:,t)=logsumexp(repmat(logalpha1(:,t-1),1,H),repmat(phatvgh1',H,1).*phghm')+1e-100;
end
loglik1 = logsumexp(logalpha1(:,T),ones(H,1)); 

if any(isnan(logalpha1(:)))
        error('logalpha1= NaN');
end;

if any(isinf(logalpha1(:)))
        error('logalpha1= Inf');
end;