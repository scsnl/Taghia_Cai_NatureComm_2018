

s = size(Lm,2);

if s>1
  temp_mcl = cat(3,Lm{:});
  temp_mcl = squeeze(temp_mcl(:,1,:));
  temp_Lcov = cat(4,Lcov{:,:,:});
  
  mean_mcl = mean(temp_mcl,2);
  nu_mcl = s./( sum(squeeze(temp_Lcov(1,1,:,:)),2) ... 
      + sum(temp_mcl.^2,2) ...
      - 2*mean_mcl.*sum(temp_mcl,2) ...
      + s*mean_mcl.^2 );
end
