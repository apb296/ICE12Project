function [resINQ,resEQ] =ResBondMarketPriceKnitro(q,Eqb,Para)
resINQ=[];
CoeffConsumptionPolicy=Eqb.CoeffConsumptionPolicy;
C=Eqb.C;
CoeffAPolicy=Eqb.CoeffAPolicy;
A=Eqb.A;
CoeffGamma=Eqb.CoeffGamma;
Gamma=Eqb.Gamma;
if Para.flagNaturalBorrowingLimit==1
phi=(-Para.S(1)/(1-q))*Para.NonZeroAdj;
else
phi=Para.AdhocBorrowingLimit;
end
count = 0;
for ctrIn=1:Para.NumIter
[CoeffConsumptionPolicy,C]=UpdateConsumptionCoeffLP(C,CoeffConsumptionPolicy,CoeffAPolicy,A,q,phi,Para);
count = count +1;
[CoeffAPolicy,A]=UpdateACoeff(CoeffConsumptionPolicy,C,A,CoeffAPolicy,q,phi,Para);
%[CoeffGamma,Gamma]=UpdateGammaCoeff(CoeffAPolicy,A,Gamma,CoeffGamma,phi,Para);
end

for ctrIn=1:Para.NumIter
[CoeffGamma,Gamma]=UpdateGammaCoeff(CoeffAPolicy,A,Gamma,CoeffGamma,phi,Para);
end
%options=optimset('Display','off','TolX',Para.ErrorTol,'MaxFunEvals',Para.NumIter*5);
%CoeffConsumptionPolicy=fsolve(@ (CoeffConsumptionPolicy) UpdateConsumptionCoeff(C,CoeffConsumptionPolicy,CoeffAPolicy,A,q,phi,Para)-CoeffConsumptionPolicy, CoeffConsumptionPolicy,options);
%CoeffAPolicy=fsolve(@ (CoeffAPolicy) UpdateACoeff(CoeffConsumptionPolicy,C,A,CoeffAPolicy,q,phi,Para)-CoeffAPolicy,CoeffAPolicy,options);
%CoeffGamma0=CoeffGamma;
 %[CoeffGamma,~,exitflag]=fsolve(@ (CoeffGamma) UpdateGammaCoeff(CoeffAPolicy,A,Gamma,CoeffGamma,phi,Para)-CoeffGamma,CoeffGamma0,options);
 %if ~(exitflag==1)
 %for ctrIn=1:Para.NumIter
%[CoeffGamma,Gamma]=UpdateGammaCoeff(CoeffAPolicy,A,Gamma,CoeffGamma,phi,Para);
% end
% end
resEQ=ResMarketClearing(CoeffGamma,Gamma, CoeffAPolicy,A,phi,Para);
Eqb.CoeffConsumptionPolicy=CoeffConsumptionPolicy;
Eqb.C=C;
Eqb.CoeffAPolicy=CoeffAPolicy;
Eqb.A=A;
Eqb.CoeffGamma=CoeffGamma;
Eqb.Gamma=Gamma;
Eqb.q=q;
Eqb.phi=phi;
end