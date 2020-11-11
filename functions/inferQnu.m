

s = size(Lm,2);

a = pa + .5*size(Lm{1},1);
for tt = 1:s
  kt = size(Lm{tt},2);
  b{tt} = pb*ones(1,kt-1) + .5*( diag(sum(Lcov{tt}(2:end,2:end,:),3))' + sum(Lm{tt}(:,2:end).^2,1) );
end


