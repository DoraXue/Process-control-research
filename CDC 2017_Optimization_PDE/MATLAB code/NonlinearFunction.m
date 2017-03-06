function dan = NonlinearFunction(A,a,gamma,betaT,Phi,X,Zstep)

TotalSign = betaT*(exp(-gamma./(1.+X)) - exp(-gamma))*Zstep;
NonlinearTerm = Phi*TotalSign;
dan = A*a + NonlinearTerm;