function [resEQ,Eqb] =ResBondMarketPrice(q,Eqb,Para)
resINQ=[];
CoeffConsumptionPolicy=Eqb.CoeffConsumptionPolicy;
C=Eqb.C;

CoeffGamma=Eqb.CoeffGamma;
Gamma=Eqb.Gamma;
if Para.flagNaturalBorrowingLimit==1
phi=(-Para.S(1)/(1-q))*Para.NonZeroAdj;
else
phi=Para.AdhocBorrowingLimit;
end

for ctrIn=1:Para.NumIter
[CoeffConsumptionPolicy,C]=UpdateConsumptionCoeff(C,CoeffConsumptionPolicy,q,phi,Para);

end

for ctrIn=1:Para.NumIter
[CoeffGamma,Gamma]=UpdateGammaCoeff(CoeffConsumptionPolicy,C,Gamma,CoeffGamma,phi,q,Para);
end

resEQ=ResMarketClearing(CoeffGamma,Gamma, CoeffConsumptionPolicy,C,phi,q,Para);
Eqb.CoeffConsumptionPolicy=CoeffConsumptionPolicy;
Eqb.C=C;

Eqb.CoeffGamma=CoeffGamma;
Eqb.Gamma=Gamma;
Eqb.q=q;
Eqb.phi=phi;
end

%CoeffAPolicy=Eqb.CoeffAPolicy;
%A=Eqb.A;
%[CoeffAPolicy,A]=UpdateACoeff(CoeffConsumptionPolicy,C,A,CoeffAPolicy,q,phi,Para);
%[CoeffGamma,Gamma]=UpdateGammaCoeff(CoeffAPolicy,A,Gamma,CoeffGamma,phi,Para);
%options=optimset('Display','off','TolX',Para.ErrorTol,'MaxFunEvals',Para.NumIter*2);
%CoeffConsumptionPolicy=fsolve(@ (CoeffConsumptionPolicy) UpdateConsumptionCoeff(C,CoeffConsumptionPolicy,CoeffAPolicy,A,q,phi,Para)-CoeffConsumptionPolicy, CoeffConsumptionPolicy,options);
%CoeffAPolicy=fsolve(@ (CoeffAPolicy) UpdateACoeff(CoeffConsumptionPolicy,C,A,CoeffAPolicy,q,phi,Para)-CoeffAPolicy,CoeffAPolicy,options);
%CoeffGamma0=CoeffGamma;
%[CoeffGamma,~,exitflag]=fsolve(@ (CoeffGamma) UpdateGammaCoeff(CoeffConsumptionPolicy,C,Gamma,CoeffGamma,phi,q,Para)-CoeffGamma,CoeffGamma0,options);
 %if ~(exitflag==1)
 %for ctrIn=1:Para.NumIter
%[CoeffGamma,Gamma]=UpdateGammaCoeff(CoeffConsumptionPolicy,C,Gamma,CoeffGamma,phi,q,Para);
% end
% end
%Eqb.CoeffAPolicy=CoeffAPolicy;
%Eqb.A=A;