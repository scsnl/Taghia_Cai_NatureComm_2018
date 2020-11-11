function VBpost = mstep_VBVAR(data,Ybar,VBprior,Gammas,barAlpha,Lambda,state)


[VBpost.mua VBpost.Lambda_a] = updata_Qa(data,Ybar,barAlpha,Gammas,Lambda,state);
[VBpost.nuN, VBpost.WN,VBpost.Lambda] = update_QLambda(data,Ybar,Gammas,VBprior,VBpost.mua,VBpost.Lambda_a,state);
[VBpost.aN,VBpost.bN,VBpost.barAlpha] = update_Qalpha(VBpost.mua,VBpost.Lambda_a,VBprior);
M = size(VBpost.Lambda,2);
digamma_args = repmat(VBpost.nuN + 1,1,M)-(1:M);
VBpost.lnDetLambda = sum(psi(0.5.*digamma_args)) + M*log(2) + log(det(VBpost.WN));


function [mua Lambda_a] = updata_Qa(data,Ybar,barAlpha,Gammas,Lambda,state)


M = size(data{1},1);
sum_mua = zeros(M^2,1);
Lambda_a = zeros(M^2);
parfor subj = 1:length(data)
    for t = 1:size(data{subj},2)
        Lambda_a = Lambda_a + Gammas{subj}(state,t) * Ybar{subj}(:,:,t)' * Lambda * Ybar{subj}(:,:,t);
        sum_mua = sum_mua + Gammas{subj}(state,t) * Ybar{subj}(:,:,t)' * Lambda * data{subj}(:,t);
    end
end
Lambda_a =  + Lambda_a + barAlpha*eye(M^2);
Ca = pinv(Lambda_a);
mua = Ca*sum_mua;



function [nuN,WN,Lambda] = update_QLambda(data,Ybar,Gammas,VBprior,mua,Lambda_a,state)

Ca = pinv(Lambda_a);
M = size(data{1},1);
Gammas_mat = cell2mat(Gammas);
Neff = sum(Gammas_mat(state,:));
nuN = Neff + VBprior.no;
invWN = zeros(M); 
parfor subj = 1:length(data)
    for t = 1:size(data{subj},2)
        Syy =   data{subj}(:,t)*data{subj}(:,t)';
        SXaX =  Ybar{subj}(:,:,t)*(Ca + mua * mua')*Ybar{subj}(:,:,t)';
        SyaX =  data{subj}(:,t) * mua' * Ybar{subj}(:,:,t)';
        SXay =  Ybar{subj}(:,:,t) * mua * data{subj}(:,t)';
        invWN = invWN + Gammas{subj}(state,t)*(Syy + SXaX - SyaX - SXay);
    end
end
invWN = pinv(VBprior.Wo) + invWN;
WN = pinv(invWN);
Lambda = nuN*WN;

function [aN,bN,barAlpha] = update_Qalpha(mua,Lambda_a,VBprior)
Ca = pinv(Lambda_a);
aN = VBprior.ao + length(mua)/2;
bN = VBprior.bo + trace(mua*mua' + Ca);
barAlpha = min(aN/bN,1e4);




