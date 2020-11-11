function  [optimalLocalDim, localDimSubj] = decideOnLocalComplexityBasedOnEnergy(data,energy)
if nargin<2
    energy = 0.9;
end
if energy>1
    error('energy is defined as: (0,1]')
end


for subj=1:length(data)
    [~,~,latent] = princomp(data{subj}');
    ix= find((cumsum(latent)/sum(latent)<energy));
    if isempty(ix)==1 
    ix= find((cumsum(latent)/sum(latent)<(energy+0.2)));
    end
    localDimSubj(subj) = ix(end);
end

optimalLocalDim = max(localDimSubj);

