function res=ResMarketClearing(CoeffGamma,Gamma, CoeffConsumptionPolicy,C,phi,q,Para)
N=Gamma(1).n+C(1).n-1;
aMin=phi;
S=Para.S;
%aMax=(Para.S(2)/(1-Para.delta))*Para.NonZeroAdj;
aMax=-phi;
%aGrid=aMin+rand(2*N,1)*(aMax-aMin);
aGrid=funnode(C(1));
%sSize=Para.sSize;
%for nn=1:N
%    H(nn,:)=aGrid(nn).^([0:N-1]);
%end
for inx_s=1:Para.sSize
ANew(:,inx_s)=min(max((aGrid+S(inx_s)-funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid))/q,aMin),aMax); % Savings given a,s
%b=(funeval(CoeffGamma(:,inx_s),Gamma(inx_s),aGrid,1)).*ANew(:,inx_s);
%Int(inx_s)=polyval(polyint(polyfit(aGrid,b,N)),-phi)-polyval(polyint(polyfit(aGrid,b,N)),phi);
Integrand=@(a) (funeval(CoeffGamma(:,inx_s),Gamma(inx_s),a,1))...
    .*min(max((a+S(inx_s)-funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),a))./q,aMin),aMax);
Int(inx_s)=quadrect(Integrand,5*N,aMin,aMax,'lege');
end
res=sum(Int)/abs(aMax);
% 
%
end
