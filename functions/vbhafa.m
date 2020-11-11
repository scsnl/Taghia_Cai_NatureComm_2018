function net = vbhafa(Y, nStates, maxdim, nIter,tol,pcaflag,approach, net)
global pa pb alpha alphaa alphapi psii u stran sprior Lm Wa Wpi Lcov...
      Xm Xcov Qns Qnss QnsCell a b logOutProbs Fhist F iter
Ycell = Y;
Y = cell2mat(Y);
[p n] = size(Y);

if nargin<2
      error('at least two inputes are expected: data, number of States')
end
if nargin<3,
      k = (p-1) + 1;
else
      k = maxdim + 1;
end
if nargin<4, nIter = 10; end
if nargin<5, tol = 1e-3; end
if nargin<6, pcaflag = 0; end
if nargin<7, approach = 1; end

if nargin==8 % given hparams, also perhaps params
      Wa = net.kmeans_init.Wa_kmeans;
      Wpi = net.kmeans_init.Wpi_kmeans;
      QnsCell = net.kmeans_init.QnsCell_kmeans;
      Lm = net.kmeans_init.Lm_init;
      net = struct( ...
            'model','variational SAFA', ...
            'hparams',struct('mcl',[],'psii',[],'pa',[],'pb',[],...
            'alpha',[],'alphaa',[],'alphapi',[]), ...
            'params',struct('Lm',[],'Lcov',[],'u',[],'stran',[],'sprior',[], 'Ybar_old', []), ...
            'hidden',struct('Xm',[],'Xcov',[],'Qns',[],'Qnss',[],'a',[],'b',[]), ...
            'kmeans_init', struct('Wa_kmeans', [], 'Wpi_kmeans', [], 'QnsCell_kmeans', []),...
            'Fhist',[], 'LL', []);
      pa = 1; pb = 1; alpha=1; alphaa = 1; alphapi = 1; psii = ones(p,1); u = 1;
      for ss=1:nStates
            Lcov{ss} = repmat(eye(k),[1 1 p]);
            mean_mcl = mean(Y,2); nu_mcl = (1./std(Y,0,2)).^2;
            Xcov{ss} = zeros(k,k,n);
      end
      s = nStates;
      pu = alpha*ones(1,s)/s;
      ua = ones(1,s)*(alphaa/s);
      upi = ones(1,s)*(alphapi/s);
      net.hidden = [];
else % nothing given, generate random everything
      net = struct( ...
            'model','variational SAFA', ...
            'hparams',struct('mcl',[],'psii',[],'pa',[],'pb',[],...
            'alpha',[],'alphaa',[],'alphapi',[]), ...
            'params',struct('Lm',[],'Lcov',[],'u',[],'stran',[],'sprior',[], 'Ybar_old', []), ...
            'hidden',struct('Xm',[],'Xcov',[],'Qns',[],'Qnss',[],'a',[],'b',[]), ...
            'kmeans_init', struct('Wa_kmeans', [], 'Wpi_kmeans', [], 'QnsCell_kmeans', []),...
            'Fhist',[], 'LL', []);
      pa = 1; pb = 1; alpha=1; alphaa = 1; alphapi = 1; psii = ones(p,1); u = 1;
      for ss=1:nStates
            Lm{ss} = randn(p,k); Lcov{ss} = repmat(eye(k),[1 1 p]);
            mean_mcl = mean(Y,2); nu_mcl = (1./std(Y,0,2)).^2;
            Xcov{ss} = zeros(k,k,n);
      end
      s = nStates;
      pu = alpha*ones(1,s)/s;
      ua = ones(1,s)*(alphaa/s);
      upi = ones(1,s)*(alphapi/s);
      F = -inf;
      Fhist = [];
      [Wa,Wpi,QnsCell] = initPoteriors(Ycell,nStates, 'subject');
      warning('on')
      Wa_kmeans = Wa;
      Wpi_kmeans = Wpi;
      QnsCell_kmeans = QnsCell;
      net.kmeans_init.Wa_kmeans = Wa_kmeans;
      net.kmeans_init.Wpi_kmeans= Wpi_kmeans;
      net.kmeans_init.QnsCell_kmeans= QnsCell_kmeans;
      net.kmeans_init.Lm_init  = Lm;
end
if nargin<9
      plotFlag=0;
end
psimin = 1e-5;


nSubjs = length(Ycell);
Qns = [];
for ns=1:nSubjs
      Qns = [Qns; QnsCell{ns}];
end
% display(['initial weights: ', mat2str(sum(Qns))]);
stran = exp(  psi(Wa) - repmat( psi(sum(Wa,2)) ,[1 s])  );
sprior = exp(  psi(Wpi) - psi(sum(Wpi,2))  );
% learning for the first iteration
iter = 1;
inferQnu;
inferQX;
inferAR3;
inferQL;
inferpsii2;
infermcl;
computeLogOutProbs;
[loglik,QnsCell,Qnss] = vbhmmEstep(Ycell,stran', sprior,logOutProbs);
Qns = [];
for ns=1:nSubjs
      Qns = [Qns; QnsCell{ns}];
end

ll = [];
F = -inf;
dF = inf;
reverseStr = '';

current_weights = -inf*sum(Qns) ;
improvement = inf;

while  iter<=nIter && improvement>tol
      learnAR_FA;
      if nargin~=8
            computeLowerBound;
            ll = [ll F+loglik];
            
            if iter==1
                  dF = Inf;
            else
                  temp = diff(ll);
                  dF =  temp(end);
            end
      end
      improvement = sum(abs(sum(Qns) - current_weights));
      current_weights = sum(Qns);
%       display(['iteration: ', num2str(iter),'.....', ' improvement:  ', num2str(improvement)]);
      msg = sprintf('iteration %d/%d', iter, nIter);
      fprintf([reverseStr, msg]);
      reverseStr = repmat(sprintf('\b'), 1, length(msg));
       iter = iter+1;
end

% normalize to account for the epsilon in the computation of HMM
sum_stran = sum(stran,2);
for state=1:nStates
      stran(state, :) = stran(state,:)/sum_stran(state);
end

net.hparams.mcl = [mean_mcl nu_mcl];
net.hparams.psii = psii;
net.hparams.pa = pa;
net.hparams.pb = pb;
net.hparams.alpha = alpha;
net.params.stran = stran;
net.params.sprior = sprior;
net.params.Wa = Wa;
net.params.Wpi = Wpi;
net.params.Lm = Lm;
net.params.Lcov = Lcov;
net.params.u = u;
net.hidden.Xm = Xm;
net.hidden.Xcov = Xcov;
net.hidden.a = a;
net.hidden.b = b;
net.hidden.Qns = Qns;
net.hidden.QnsCell = QnsCell;
net.hidden.Qnss = Qnss;
net.Fhist = Fhist;
net.LL= ll;
net.logOutProbs = logOutProbs;
net.ARmodel = ARpost;
net.params.Ybar_old = Ybar_old;
net.hparams.pcaflag = pcaflag;


