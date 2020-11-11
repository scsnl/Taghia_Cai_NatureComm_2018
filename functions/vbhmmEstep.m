function [loglik,phtgV1T,phthtpgV1T] = vbhmmEstep(v,stran,sprior,logOutputs)

nSubjs = length(v);
parfor ns=1:nSubjs
      [logalpha,loglik1]=VBHMMforward(v{ns},stran,sprior,logOutputs{ns});
      logbeta=VBHMMbackward(v{ns},stran,logOutputs{ns});
      [phtgV1T1,phthtpgV1T1]=VBHMMsmooth(logalpha,logbeta,logOutputs{ns},stran,v{ns}); 
      for t = 1:size(phthtpgV1T1,3)
            phthtpgV1T1(:,:,t) = phthtpgV1T1(:,:,t)';
      end
      loglik{ns} = loglik1;
      phtgV1T{ns} = phtgV1T1'; 
      phthtpgV1T{ns} = phthtpgV1T1;  
end
loglik = sum(cell2mat(loglik));
