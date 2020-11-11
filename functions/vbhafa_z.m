function net = vbhafa_z(Y, Qns, maxdim, nIter, tol, pcaflag)
QnsCell = {Qns};
Ycell = Y;
Y = cell2mat(Y);
nStates = size(Qns, 2);
[p n] = size(Y);
k = maxdim + 1;
approach = 1;
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
F = -inf;
Fhist = [];
nSubjs = length(Ycell);
% learning for the first iteration
iter = 1;
psimin = 1e-5;
inferQnu;
inferQX;
inferAR3;
inferQL;
inferpsii2;
infermcl;
computeLogOutProbs;
ll = [];
F = -inf;
reverseStr = '';
current_weights = -inf*sum(Qns) ;
improvement = inf;
while  iter<=nIter && improvement>tol
      learnAR_FA_z;
      computeLowerBound_z;
      ll = [ll F];
      msg = sprintf('iteration %d/%d', iter, nIter);
      fprintf([reverseStr, msg]);
      reverseStr = repmat(sprintf('\b'), 1, length(msg));
       iter = iter+1;
end
% getAREstimates;

net.hparams.mcl = [mean_mcl nu_mcl];
net.hparams.psii = psii;
net.hparams.alpha = alpha;
net.params.Lm = Lm;
net.params.Lcov = Lcov;
net.params.u = u;
net.hidden.Xm = Xm;
net.hidden.Xcov = Xcov;
net.hidden.a = a;
net.hidden.b = b;
net.Fhist = Fhist;
net.LL= ll;
net.logOutProbs = logOutProbs;
net.ARmodel = ARpost;
net.params.Ybar_old = Ybar_old;
net.hidden.QnsCell = QnsCell;


