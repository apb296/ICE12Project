function [CoeffConsumptionPolicy,C]=UpdateConsumptionCoeff(C,CoeffConsumptionPolicy,q,phi,Para)
sigma=Para.sigma;
delta=Para.delta;
S=Para.S;
P=Para.P;
aMin=phi;
%aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aMax=-phi;
sSize=Para.sSize;
GridDensity=Para.GridDensity;
OrderOfApproxConsumptionPolicy=Para.OrderOfApproxConsumptionPolicy;
ApproxMethod=Para.ApproxMethod;
aGridSize=GridDensity*OrderOfApproxConsumptionPolicy;
aGrid=funnode(C(1));
CNew=ones(aGridSize,sSize);
ANew=ones(aGridSize,sSize);
% ANew = A(a,s) given q. Savings given state today
CTomorrow=ones(aGridSize,sSize);
% CTomorrow(s'| a,s) = C[A(a,s),s']
    
for inx_s=1:sSize % state today - s
     
    if ~(C(inx_s).a==phi)
        C(inx_s) = fundefn(ApproxMethod,OrderOfApproxConsumptionPolicy ,aMin,aMax);
    end
    %MaxConsumptionToday(:,inx_s)=aGrid+repmat(S(inx_s),aGridSize,1)-repmat(q*aMin,aGridSize,1)*(1/Para.NonZeroAdj);
    %ConsumptionToday(:,inx_s)=min(funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid),MaxConsumptionToday(:,inx_s));
    ConsumptionToday(:,inx_s)=funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid);
    ANew(:,inx_s)=(aGrid+repmat(S(inx_s),aGridSize,1)-ConsumptionToday(:,inx_s))./q;
    ANew(:,inx_s)=max(min(ANew(:,inx_s),aMax),aMin);
    for inx_sTomorrow=1:sSize % state tomorrow - s'
        CTomorrow(:,inx_sTomorrow)=max(funeval(CoeffConsumptionPolicy(:,inx_sTomorrow),C(inx_sTomorrow),ANew(:,inx_s)),.0001); % consumption tomorrow (s'|a,s)
    end
    MuTomorrow=CTomorrow.^(-sigma);
    EMuTomorrow=MuTomorrow*P(inx_s,:)';
    MuToday(:,inx_s)=delta*EMuTomorrow./q;
    %CNew(:,inx_s)=min(MuToday(:,inx_s).^(-1/sigma),MaxConsumptionToday(:,inx_s));
    CNew(:,inx_s)=MuToday(:,inx_s).^(-1/sigma);
    CoeffConsumptionPolicy(:,inx_s)=funfitxy(C(inx_s),aGrid,CNew(:,inx_s));
end






