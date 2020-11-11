

n = size(Y,2);
p = size(Y,1);

for tt = 1:s
  kt = size(Lm{tt},2);
  T1 = reshape(reshape(Lcov{tt}(2:end,2:end,:),(kt-1)*(kt-1),p)*psii,kt-1,kt-1) ...
      + Lm{tt}(:,2:end)'*diag(psii+eps)*Lm{tt}(:,2:end);
  Xcov{tt}              = zeros(kt,kt);
  Xcov{tt}(2:end,2:end) = inv( eye(kt-1)+T1 );
  trXm{tt}              = Xcov{tt}(2:end,2:end)*Lm{tt}(:,2:end)'*diag(psii+eps)*(Y-Lm{tt}(:,1)*ones(1,n));
  Xm{tt}                = [ones(1,n); trXm{tt}];
end


