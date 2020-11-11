function [counter,sprior,posteriors] = initPoteriors(data, nStates, method_kmeans)
warning('off')
if nargin<3
      method_kmeans = 'subject';
end
switch method_kmeans
      case 'group'
            nSubjs = length(data);
            [M,T] = size(cell2mat(data));
            states = kmeans(cell2mat(data)',nStates,'Replicates',10,'Distance','sqEuclidean','emptyaction','drop');
            for hh = 1:nStates
                  ix = states == hh;
                  temp(hh,:) = ix;
            end
            posteriors = temp'; 
            
            sprior = zeros(1,nStates);
            for hh = 1:nStates
                  ix = states == hh;
                  sprior(hh) = sum(ix);
            end
            counter = zeros(nStates);
            for t1 = 1:T-1
                  t2 = t1+1;
                  m = states(t2); n = states(t1);
                  counter(m,n) = counter(m,n) + 1;
            end
            counter = counter';
      case 'subject'
            
            nSubjs = length(data);
            [M,T] = size(cell2mat(data));
            parfor ns = 1:nSubjs
                  statesCell{ns} = kmeans(data{ns}',nStates,'Replicates',10,'Distance','sqEuclidean','emptyaction','drop');
            end
            states = cell2mat(statesCell);
            states = states(:);
            for ns = 1:nSubjs
                  for hh = 1:nStates
                        ix = statesCell{ns} == hh;
                        temp(hh,:) = ix;
                  end
                  posteriors{ns} = temp';
            end
            
            sprior = zeros(1,nStates);
            for hh = 1:nStates
                  ix = states == hh;
                  sprior(hh) = sum(ix);
            end
            counter = zeros(nStates);
            for t1 = 1:T-1
                  t2 = t1+1;
                  m = states(t2); n = states(t1);
                  counter(m,n) = counter(m,n) + 1;
            end
            counter = counter';
end