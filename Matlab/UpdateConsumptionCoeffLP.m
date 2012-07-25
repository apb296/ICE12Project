function [C,CoeffConsumptionPolicy]=UpdateConsumptionCoeffLP(C,CoeffConsumptionPolicy,CoeffAPolicy,A,q,phi,Para);
% Coefficient update by LP with shape constraints
%   Detailed explanation goes here

sigma=Para.sigma;
delta=Para.delta;
S=Para.S;
P=Para.P;
aMin=phi;
aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
sSize=Para.sSize;
GridDensity=Para.GridDensity;
OrderOfApproxConsumptionPolicy=Para.OrderOfApproxConsumptionPolicy;
ApproxMethod='cheb';
%aGridSize=GridDesity*OrderOfApproxConsumptionPolicy;

aGrid=funnode(C(1));%using Chebychef nodes for fitting
aGridSize=length(aGrid);

CNew=ones(aGridSize,sSize);
ANew=ones(aGridSize,sSize);
% ANew = A(a,s) given q. Savings given state today
CTomorrow=ones(aGridSize,sSize);
% CTomorrow(s'| a,s) = C[A(a,s),s']
    
for inx_s=1:sSize % state today - s
    
    ANew(:,inx_s)=max(min(funeval(CoeffAPolicy(:,inx_s),A(inx_s),aGrid),aMax),aMin); % Savings given a,s
    
    if ~(C(inx_s).a==phi)
        C(inx_s) = fundefn(ApproxMethod,OrderOfApproxConsumptionPolicy ,aMin,aMax);
    end
    
    for inx_sTomorrow=1:sSize % state tomorrow - s'
        CTomorrow(:,inx_sTomorrow)=max(funeval(CoeffConsumptionPolicy(:,inx_sTomorrow),C(inx_sTomorrow),ANew(:,inx_s)),.001); % consumption tomorrow (s'|a,s)
    end
    CNew(:,inx_s)=(delta*CTomorrow.^(-sigma)*P(inx_s,:)'/q).^(-1/sigma);
    CoeffConsumptionPolicy(:,inx_s)=funfitxy(C(inx_s),aGrid,CNew(:,inx_s));

end


end

