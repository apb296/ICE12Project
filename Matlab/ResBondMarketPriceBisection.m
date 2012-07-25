function [resEQ] =ResBondMarketPriceBisection(q,Eqb,Para,OutputFlag)
%q
%global Eqb;
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

ctrIn=1;
DiffCoeffError=1;
CoeffConsumptionPolicy00=CoeffConsumptionPolicy;
while ctrIn<Para.NumIter && DiffCoeffError>Para.ErrorTol*10
CoeffConsumptionPolicyOld=CoeffConsumptionPolicy;
[CoeffConsumptionPolicy,C]=UpdateConsumptionCoeff(C,CoeffConsumptionPolicy,q,phi,Para);
DiffCoeffError=sum(sum(abs(CoeffConsumptionPolicy-CoeffConsumptionPolicyOld)));
ctrIn=ctrIn+1;
end


 if ctrIn==Para.NumIter
     DiffCoeffError;
 options=optimset('Display','off','TolX',Para.ErrorTol*10,'MaxFunEvals',Para.NumIter*2);
  [CoeffConsumptionPolicy,~,exitflag]=fsolve(@ (CoeffConsumptionPolicy) UpdateConsumptionCoeff(C,CoeffConsumptionPolicy,q,phi,Para)-CoeffConsumptionPolicy, CoeffConsumptionPolicy00,options);
 % if ~(exitflag==1)
 %    CoeffConsumptionPolicy=CoeffConsumptionPolicy00;
 %end
 disp('using newton for consumption policy')
 end




ctrIn=1;
DiffGammaError=1;
while ctrIn<Para.NumIter && DiffGammaError>Para.ErrorTol*10
CoeffGammaOld=CoeffGamma;
[CoeffGamma,Gamma]=UpdateGammaCoeff(CoeffConsumptionPolicy,C,Gamma,CoeffGamma,phi,q,Para);
DiffGammaError=sum(sum(abs(CoeffGamma-CoeffGammaOld)));
ctrIn=ctrIn+1;
end

resEQ=ResMarketClearing(CoeffGamma,Gamma, CoeffConsumptionPolicy,C,phi,q,Para);
Eqb.CoeffConsumptionPolicy=CoeffConsumptionPolicy;
Eqb.C=C;
Eqb.CoeffGamma=CoeffGamma;
Eqb.Gamma=Gamma;
Eqb.q=q;
Eqb.phi=phi;

% if strcmpi(OutputFlag,'solver')==1
%     varargout{1}=[];
% else
%     varargout{1}=Eqb;
% end

end

%CoeffAPolicy=Eqb.CoeffAPolicy;
%A=Eqb.A;
%[CoeffAPolicy,A]=UpdateACoeff(CoeffConsumptionPolicy,C,A,CoeffAPolicy,q,phi,Para);
%[CoeffGamma,Gamma]=UpdateGammaCoeff(CoeffAPolicy,A,Gamma,CoeffGamma,phi,Para);
%options=optimset('Display','off','TolX',Para.ErrorTol,'MaxFunEvals',Para.NumIter*5);
%CoeffConsumptionPolicy=fsolve(@ (CoeffConsumptionPolicy) UpdateConsumptionCoeff(C,CoeffConsumptionPolicy,CoeffAPolicy,A,q,phi,Para)-CoeffConsumptionPolicy, CoeffConsumptionPolicy,options);
%CoeffAPolicy=fsolve(@ (CoeffAPolicy) UpdateACoeff(CoeffConsumptionPolicy,C,A,CoeffAPolicy,q,phi,Para)-CoeffAPolicy,CoeffAPolicy,options);
%CoeffGamma0=CoeffGamma;
% [CoeffGamma,~,exitflag]=fsolve(@ (CoeffGamma) UpdateGammaCoeff(CoeffConsumptionPolicy,C,Gamma,CoeffGamma,phi,q,Para)-CoeffGamma,CoeffGamma0,options);
% if ~(exitflag==1)
% for ctrIn=1:Para.NumIter
%[CoeffGamma,Gamma]=UpdateGammaCoeff(CoeffConsumptionPolicy,C,Gamma,CoeffGamma,phi,q,Para);
% end
% end
%Eqb.CoeffAPolicy=CoeffAPolicy;
%Eqb.A=A;