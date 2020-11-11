function stats = posthocVARfromData(data, net, estStates, sortFlag, nIter, p_ref)
if nargin<5
      nIter =100;
      p_ref = 0.01;
end
if nargin<6
      p_ref = 0.01;
end

K = length( net.ARmodel);

M = size(data{1},1);
VBprior.ao = 10^-3; VBprior.bo = 10^-3;
VBprior.no = M + 2;
VBprior.Wo = 10^3*eye(M);

noSubjs = length(data);
Ybar = cell(1,noSubjs);
parfor subj = 1:noSubjs
      Ybar{subj}(:,:,1) = zeros(M,M^2);
      for t = 1:size(data{subj},2)-1
            Ybart = zeros(M,M^2);
            st = 1;
            for m = 1:M
                  yt_1 = data{subj}(:,t);
                  Ybart(m,st:st+M-1) = yt_1';
                  st = st + M;
            end
            Ybar{subj}(:,:,t+1) = Ybart;
      end
end
for subj = 1:length(data)
      Gammas{subj} = net.hidden.QnsCell{subj}';
end
reverseStr = '';
for iter = 1:nIter
      
      for k = 1:K
            if iter == 1
                  barAlpha = (VBprior.ao)/(VBprior.bo); 
                  Lambda = VBprior.Wo;
            else
                  barAlpha = VBpost{k}.barAlpha;
                  Lambda = VBpost{k}.Lambda;
            end
            VBpost{k} = mstep_VBVAR(data,Ybar,VBprior,Gammas,barAlpha,Lambda,k);
      end
      msg = sprintf('iteration %d/%d', iter, nIter);
      fprintf([reverseStr, msg]);
      reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
display(' ')
if sortFlag==1
      uniqueStates = unique(estStates);
      for i=1:numel(uniqueStates)
            countStates(i) = sum(estStates==uniqueStates(i));
      end
[~, id]   = sort(countStates,'descend');
dominant_states = uniqueStates(id);
elseif sortFlag==0
      dominant_states = unique(estStates);
else
      error('sortFlag either on or off')
end

pval_bon = p_ref/(length(dominant_states )*M*(M-1));


thresh = norminv(1-pval_bon);
Map = zeros(M,M,length(dominant_states)); Theta_A = Map; Pvals = Map;
for k = 1:numel(dominant_states)
      A = reshape(VBpost{dominant_states(k)}.mua,M,M)';
      Ca = pinv(VBpost{dominant_states(k)}.Lambda_a);
      diag_Ca = reshape(diag(Ca)',M,M);
      Theta_A(:,:,k) = A./sqrt(diag_Ca);
      Map(:,:,k) = double(abs(Theta_A(:,:,k)) >= thresh);
      Pvals(:,:,k) = 1-normcdf(abs(Theta_A(:,:,k)));
end

stats.Pvals = Pvals;
stats.Acoeffs = A;
stats.Zscores =  Theta_A;
stats.Map = Map;

