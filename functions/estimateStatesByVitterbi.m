function [states, stateCell] =  estimateStatesByVitterbi(data,model,logOutProbs)

Wa = model.Wa;
Wa = Wa';
K = size(Wa,1);
Aest = Wa./repmat(sum(Wa,2),1,K); 
Aest = Aest';
Wpi = model.Wpi;
piest = Wpi./sum(Wpi);
nSubjs = length(data);
states = [];
for ns = 1:nSubjs
    [M,T] = size(data{ns});
    pvgh = zeros(K,T);
    for h = 1:K
        for t = 1:1
            pvgh(h,t) =  exp(logOutProbs{ns}(h,t));
        end
        for t = 2:T
            pvgh(h,t) =  exp(logOutProbs{ns}(h,t));
        end
    end
    states = [states viterbi_path(piest, Aest',pvgh)];
    stateCell{ns} = viterbi_path(piest, Aest',pvgh);
end

