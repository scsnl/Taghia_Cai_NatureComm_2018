
n = size(Y,2);
p = size(Y,1);
s = size(Lm,2);

for tt = 1:s
  kt = size(Lm{tt},2);
  mean_Lambda = [mean_mcl zeros(p,kt-1)];
  num{tt} = [nu_mcl repmat(a./b{tt},[p 1]) ];
  temp = Xm{tt}.*repmat(Qns(:,tt)',kt,1);
  T2 = Xcov{tt}*sum(Qns(:,tt),1) + Xm{tt}*temp';
  T3 = diag(psii)* Y *temp';
  for q = 1:p
    Lcov{tt}(:,:,q) = inv( diag(num{tt}(q,:)) + psii(q)*T2);
    Lm{tt}(q,:) = ( T3(q,:) + mean_Lambda(q,:)*diag(num{tt}(q,:)) ) * Lcov{tt}(:,:,q);
  end
end


