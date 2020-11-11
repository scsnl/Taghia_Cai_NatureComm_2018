n = size(Y,2);
p = size(Y,1);
s = nStates;
ldim = size(Xm{1},1)-1;
for subj = 1:nSubjs
      Gammas{subj} = QnsCell{subj}';
end

for tt=1:s
      parfor ns=1:nSubjs
            XmState{ns}(:,:) = preprocess(Xm{tt}(2:end,n*(ns-1)/nSubjs+1:(n*ns/nSubjs)));
      end
      Ybar = data_arrange(XmState);
      if sum(sum(any(isnan(cell2mat(Ybar)))))
            Ybar = Ybar_old;
      end
      if sum(sum(any(isinf(cell2mat(Ybar)))))
            Ybar = Ybar_old;
      end
      if sum(sum(any(isnan(cell2mat(XmState)))))
            XmState = XmState_old;
      end
      if sum(sum(any(isinf(cell2mat(XmState)))))
            XmState = XmState_old;
      end
      ARprior = set_ARhyperpriors(XmState);
      
      if iter ==1
            parfor tt2=1:nStates
                  barAlpha = (ARprior.ao)/(ARprior.bo);
                  Lambda = ARprior.Wo;
                  ARpost{tt2} = mstep_VBVAR(XmState,Ybar,ARprior,Gammas,barAlpha,Lambda,tt2);
            end
      else
            parfor tt2 = 1:nStates
                  barAlpha = ARpost{tt2}.barAlpha;
                  Lambda = ARpost{tt2}.Lambda;
                  ARpost{tt2} = mstep_VBVAR(XmState,Ybar,ARprior,Gammas,barAlpha,Lambda,tt2);
            end
      end
      B_bar = reshape(ARpost{tt}.mua,ldim,ldim)'; 
      CovNoise = ARpost{tt}.nuN.*ARpost{tt}.WN;
      Ex{tt}(:,1) = Xm{tt}(2:end,1);
      CovXi{tt}(:,:,1) = Xcov{tt}(2:end,2:end);
      %...........................................
      switch approach
            case 1
                  for nn=2:n
                        Ex{tt}(:,nn) = B_bar * Xm{tt}(2:end,nn-1);
                        EXXj = Xcov{tt}(2:end,2:end) + Xm{tt}(2:end,nn-1)*Xm{tt}(2:end,nn-1)';
                        CovBXj = B_bar* EXXj* B_bar' - Ex{tt}(:,nn)*Ex{tt}(:,nn)';
                        CovXi{tt}(:,:,nn) = CovBXj  + CovNoise;
                        [~,pchol] = chol(CovXi{tt}(:,:,nn));
                        if pchol~=0
                              CovXi{tt}(:,:,nn) = CovXi{tt}(:,:,nn-1);
                        end
                  end
                  Xm{tt}(2:end,:) = Ex{tt};
                  Xcov{tt}(2:end,2:end) = mean(CovXi{tt},3);
                  
            case 2
                  
                  for nn=2:n
                        Ex{tt}(:,nn) = B_bar * Xm{tt}(2:end,nn-1);
                        CovXi{tt}(:,:,nn) = B_bar* Xcov{tt}(2:end,2:end)* B_bar' + CovNoise;
                        [~,pchol] = chol(CovXi{tt}(:,:,nn));
                        if pchol~=0
                              CovXi{tt}(:,:,nn) = CovXi{tt}(:,:,nn-1);
                        end
                  end
                  Xm{tt}(2:end,:) = Ex{tt};
                  Xcov{tt}(2:end,2:end) = mean(CovXi{tt},3);
            case 3
                  
                  for nn=2:n
                        Ex{tt}(:,nn) = B_bar * Xm{tt}(2:end,nn-1);
                        CovXi{tt}(:,:,nn) = B_bar* Xcov{tt}(2:end,2:end)* B_bar' + CovNoise;
                        [~,pchol] = chol(CovXi{tt}(:,:,nn));
                        if pchol~=0
                              CovXi{tt}(:,:,nn) = CovXi{tt}(:,:,nn-1);
                        end
                  end
                  Ex{tt}(:,1) = Ex{tt}(:,2);
                  CovXi{tt}(:,:,1) = CovXi{tt}(:,:,2);
                  Xm{tt}(2:end,:) = Ex{tt};
                  Xcov{tt}(2:end,2:end) = mean(CovXi{tt},3);
                  
            otherwise
                  error('choose a value between 1-3 for approach')
      end
      
end
