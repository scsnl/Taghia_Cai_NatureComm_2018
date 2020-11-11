

n = size(Y,2);
p = size(Y,1);
s = size(Lm,2);

psi2 = zeros(p); zeta = 0;

for tt = 1:s
  kt = size(Lm{tt},2);
  temp_alt = Xm{tt}.*repmat(Qns(:,tt)',kt,1);
  temp = Xcov{tt}*sum(Qns(:,tt),1)  +  Xm{tt}*temp_alt';
  if pcaflag == 0
    psi2 = psi2 + (repmat(Qns(:,tt)',p,1).*Y) * (Y-2*Lm{tt}*Xm{tt})' ...
	+ Lm{tt}*temp*Lm{tt}';
    for q = 1:p
      psi2(q,q) = psi2(q,q) + trace( Lcov{tt}(:,:,q)*temp );
    end
  else
    zeta = zeta + ( + sum((Y.*(Y-2*Lm{tt}*Xm{tt}))*Qns(:,tt),1) ...
	+ trace(Lm{tt}*temp*Lm{tt}') ...
	+ trace(sum(Lcov{tt},3)*temp) ...
	);
  end
end



if pcaflag == 0
  psi2 = 1/n * psi2;
  psii = 1./diag(psi2);
  if any( (1./psii) < psimin )
    fprintf('psi threshold');
    psii(find(psii>(1./psimin))) = 1./psimin;
  end
else
  psii = n*p*ones(p,1)*zeta^-1;
end

