nStates = size(Lm,2);
alphaa = 1;
alphapi = 1;
ua = ones(1,s)*(alphaa/s);
upi = ones(1,s)*(alphapi/s);

FPI=0;
for state = 1:nStates,
    FPI = FPI - kldirichlet(Wa(state,:),ua);
end;
Fpi = - kldirichlet(Wpi,upi);
lowerBoundHMMparams = FPI + Fpi;