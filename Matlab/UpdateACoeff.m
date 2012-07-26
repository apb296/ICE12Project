function [CoeffAPolicy,A]=UpdateACoeff(CoeffConsumptionPolicy,C,A,CoeffAPolicy,q,phi,Para)
sigma=Para.sigma;
delta=Para.delta;
S=Para.S;
P=Para.P;
aMin=phi;
aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
sSize=Para.sSize;
GridDensity=Para.GridDensity;
OrderOfApproxAPolicy=Para.OrderOfApproxAPolicy;
ApproxMethod='cheb';
aGridSize=GridDensity*OrderOfApproxAPolicy;
aGrid=funnode(A(1));
ANew=ones(aGridSize,sSize);
for inx_s=1:sSize
ANew(:,inx_s)=max(min((aGrid+repmat(S(inx_s),aGridSize,1)-funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid))./q,aMax),aMin);
end
for inx_s=1:sSize
    
    if ~(A(inx_s).a==phi)
        A(inx_s) = fundefn(ApproxMethod,OrderOfApproxAPolicy ,aMin,aMax);    

    end
    

CoeffAPolicy(:,inx_s)=funfitxy(A(inx_s),aGrid,ANew(:,inx_s));
end

end

