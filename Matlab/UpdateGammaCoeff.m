function [CoeffGamma,Gamma]=UpdateGammaCoeff(CoeffConsumptionPolicy,C,Gamma,CoeffGammaOld,phi,q,Para)
sigma=Para.sigma;
delta=Para.delta;
S=Para.S;
P=Para.P;
aMin=phi;
%aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aMax=-phi;
sSize=Para.sSize;
GridDensity=Para.GridDensity;
OrderOfApproxGamma=Para.OrderOfApproxGamma;
ApproxMethod='cheb';
aGridSize=GridDensity*OrderOfApproxGamma*2;
aGrid=linspace(aMin,aMax,aGridSize)';
CoeffGamma=ones(OrderOfApproxGamma,sSize);
GammaNew=ones(aGridSize,sSize);
OrderOfApproxAInv=Para.OrderOfApproxConsumptionPolicy;
for inx_s=1:sSize
    if ~(Gamma(inx_s).a==phi)
   Gamma(inx_s) = fundefn(ApproxMethod,OrderOfApproxGamma ,aMin,aMax);
    end
   ANew(:,inx_s)=min(max((aGrid+S(inx_s)-funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid))/q,aMin),aMax); % Savings given a,s
   AInv(inx_s) = fundefn(ApproxMethod,OrderOfApproxAInv,aMin,aMax);
   CoeffAInv(:,inx_s)=funfitxy(AInv(inx_s),ANew(:,inx_s),aGrid); 
   Ainv(:,inx_s)=min(max(funeval(CoeffAInv(:,inx_s),AInv(inx_s),aGrid),aMin),aMax);
   GammaNew(:,inx_s)=max(min(funeval(CoeffGammaOld(:,inx_s),Gamma(inx_s),Ainv(:,inx_s)),1),0);
   
end
for inx_s=1:sSize
    GammaFit(:,inx_s)=GammaNew*P(inx_s,:)';
    CoeffGamma(:,inx_s)=FitGammaCoeff(GammaFit(:,inx_s),aGrid,Gamma(inx_s),phi,Para);
    %CoeffGamma(:,inx_s)=FitGammaCoeffLP(GammaFit(:,inx_s),aGrid,Gamma(inx_s),phi,Para);
end
end

