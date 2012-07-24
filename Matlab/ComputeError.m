function [Error]=ComputeError(Eqb,Para)
CoeffConsumptionPolicy=Eqb.CoeffConsumptionPolicy;
C=Eqb.C;
CoeffAPolicy=Eqb.CoeffAPolicy;
A=Eqb.A;
CoeffGamma=Eqb.CoeffGamma;
Gamma=Eqb.Gamma;
phi=Eqb.phi;
q=Eqb.q;
aMin=phi;
aMax=(Para.S(2)/(1-Para.delta))*Para.NonZeroAdj;
Error.ConsumptionPolicy=0;
Error.APolicy=0;
Error.Gamma=0;
Error.MarketClearing=0;
[CoeffConsumptionPolicyNew,C]=UpdateConsumptionCoeff(C,CoeffConsumptionPolicy,CoeffAPolicy,A,q,phi,Para);
[CoeffAPolicyNew,A]=UpdateACoeff(CoeffConsumptionPolicy,C,A,CoeffAPolicy,q,phi,Para);
[CoeffGammaNew,Gamma]=UpdateGammaCoeff(CoeffAPolicy,A,Gamma,CoeffGamma,phi,Para);
res=ResMarketClearing(CoeffGamma,Gamma, CoeffAPolicy,A,phi,Para);
NumTestPoints=10;
TestaGrid=aMin+rand(NumTestPoints,1)*(aMax-aMin);
for inx_s=1:Para.sSize
 Error.ConsumptionPolicy= Error.ConsumptionPolicy./(NumTestPoints)+sqrt(sum((funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),TestaGrid)-funeval(CoeffConsumptionPolicyNew(:,inx_s),C(inx_s),TestaGrid)).^2 ) )/(NumTestPoints);
 Error.APolicy= Error.APolicy/(NumTestPoints)+sqrt(sum((funeval(CoeffAPolicy(:,inx_s),A(inx_s),TestaGrid)-funeval(CoeffAPolicyNew(:,inx_s),A(inx_s),TestaGrid)).^2 ) )/(NumTestPoints);
 Error.Gamma= Error.Gamma/(NumTestPoints)+sqrt(sum((funeval(CoeffGamma(:,inx_s),Gamma(inx_s),TestaGrid)-funeval(CoeffGammaNew(:,inx_s),Gamma(inx_s),TestaGrid)).^2 ) )/(NumTestPoints);
 Error.MarketClearing=abs(res);
end
end