function res=ResMarketClearing(CoeffGamma,Gamma, CoeffAPolicy,A,phi,Para)
N=Gamma(1).n+A(1).n-1;
aMin=phi;
aMax=(Para.S(2)/(1-Para.delta))*Para.NonZeroAdj;
%aGrid=aMin+rand(2*N,1)*(aMax-aMin);
%aGrid=linspace(aMin,aMax,3*N)';
%sSize=Para.sSize;
%for nn=1:N
%    H(nn,:)=aGrid(nn).^([0:N-1]);
%end
for inx_s=1:Para.sSize
 %b=(funeval(CoeffGamma(:,inx_s),Gamma(inx_s),aGrid,1)).*(funeval(CoeffAPolicy(:,inx_s),A(inx_s),aGrid));
%Int(inx_s)=polyval(polyint(polyfit(aGrid,b,N)),-phi)-polyval(polyint(polyfit(aGrid,b,N)),phi);
Integrand=@(a) (funeval(CoeffGamma(:,inx_s),Gamma(inx_s),a,1)).*(funeval(CoeffAPolicy(:,inx_s),A(inx_s),a));
Int(inx_s)=quadrect(Integrand,3*N,aMin,aMax,'lege');
end
res=sum(Int);
% 
%
end
