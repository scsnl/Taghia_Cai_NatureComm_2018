function [phtgV1T,phthtpgV1T1]=VBHMMsmooth(logalpha,logbeta,logOutput,phghm,v)

T = size(v,2);  
H = size(logOutput,1);  
for t=1:T
    logphtgV1T(:,t)=logalpha(:,t)+logbeta(:,t); 
    phtgV1T(:,t)=condexp(logphtgV1T(:,t));
end
for t=2:T
    atmp=condexp(logalpha(:,t-1));
    btmp=condexp(logbeta(:,t));
    phatvgh1 = zeros(H,1);
    for h = 1:H
        phatvgh1(h) = exp(logOutput(h,t))+1e-100; % eps
    end
    phatvgh1=condexp(phatvgh1);
    phghmt=phghm;
    ctmp1 = repmat(atmp,1,H).*phghmt'.*repmat(phatvgh1'.*btmp',H,1) + 1e-100;
    phthtpgV1T1(:,:,t-1)=ctmp1./sum(sum(ctmp1));
end

if any(isnan(phthtpgV1T1(:)))
        error('phthtpgV1T1= NaN');
end;

if any(isinf(phthtpgV1T1(:)))
        error('phthtpgV1T1= Inf');
end;

if any(isnan(phtgV1T(:)))
        error('phtgV1T= NaN');
end;

if any(isinf(phtgV1T(:)))
        error('phtgV1T= Inf');
end;
