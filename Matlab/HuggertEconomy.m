% This file solves the Huggert model with projection methods 
clc
clear all
close all
%% SET PARAMETERS 
% Technology and Prefernces
sigma=5; % Risk Aversion
delta=.9; % Time discount factor
S=[.1 .2 ]; % Endowment shocks
sSize=length(S);
alpha=.5; % Probability of staying in the same state
P=[alpha 1-alpha;1-alpha alpha]; % Stochastic shock matrix
Para.sigma=sigma;
Para.delta=delta;
Para.S=S;
Para.sSize=sSize;
Para.P=P;
Para.flagNaturalBorrowingLimit=1;
Para.AdhocBorrowingLimit=(-S(1)/(1-delta))*.1;

% others 
ApproxMethod='cheb';
OrderOfApproxConsumptionPolicy=5;
OrderOfApproxAPolicy=5;
OrderOfApproxGamma=20;
ErrorTol=1e-5;
GridDensity=1;
NonZeroAdj=.95;
Para.ApproxMethod=ApproxMethod;
Para.OrderOfApproxConsumptionPolicy=OrderOfApproxConsumptionPolicy;
Para.OrderOfApproxAPolicy=OrderOfApproxAPolicy;
Para.GridDensity=GridDensity;
Para.OrderOfApproxGamma=OrderOfApproxGamma;
Para.NonZeroAdj=NonZeroAdj;
grelax=.9;
NumIter=40;
Para.NumIter=NumIter;
Para.ErrorTol=ErrorTol;
% Set Paths
if strcmp(computer,'PCWIN')
    sl='\';
else
    sl='/';
end
compeconpath=[pwd '/lib/compecon2011' sl];
texpath= [pwd sl 'Tex' sl] ;
plotpath= [pwd sl 'Graphs' sl] ;
datapath=[pwd sl 'Data' sl] ;
mkdir(texpath)
mkdir(plotpath)
mkdir(datapath)
addpath(genpath(compeconpath))
%% INITIALIZE EQB OBJECTS
% We are solving for 5 equilibrium objects 
% 1) q - price of the bond
% 2) phi - Natural borrowing limit
% 3) C[a,s,q] - Consumption Rule given q
% 4) A[a,s,q] - Savings Rule given q
% 5) Gamma[a,s] - Ergodic distribution over a,s

% -------
q=delta; % Price of the bond with complete markets
% -------


%--phi--
if Para.flagNaturalBorrowingLimit==1
phi=(-Para.S(1)/(1-q))*Para.NonZeroAdj;
else
phi=Para.AdhocBorrowingLimit;
end


% -- C[a,s,q] --
aMin=phi;
aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aGridSize=OrderOfApproxConsumptionPolicy;
C(1) = fundefn(ApproxMethod,OrderOfApproxConsumptionPolicy ,aMin,aMax);
aGrid=funnode(C(1));
C0=repmat((1-q).*aGrid,1,sSize)+repmat(S,aGridSize,1);
CoeffConsumptionPolicy=ones(OrderOfApproxConsumptionPolicy,sSize);
for inx_s=2:sSize
    C(inx_s) = fundefn(ApproxMethod,OrderOfApproxConsumptionPolicy ,aMin,aMax);

    CoeffConsumptionPolicy(:,inx_s)=funfitxy(C(inx_s),aGrid,C0(:,inx_s));
end

% -- A[a,s,q] --
aMin=phi;
aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aGridSize=GridDensity*OrderOfApproxAPolicy;
A(1) = fundefn(ApproxMethod,OrderOfApproxAPolicy ,aMin,aMax);
aGrid=funnode(A(1));
A0=repmat(aGrid,1,sSize);
CoeffAPolicy=ones(OrderOfApproxAPolicy,sSize);
for inx_s=2:sSize
    A(inx_s) = fundefn(ApproxMethod,OrderOfApproxAPolicy ,aMin,aMax);
    
    CoeffAPolicy(:,inx_s)=funfitxy(A(inx_s),aGrid,A0(:,inx_s));
end

% -- Gamma[a,s,q]---
aMin=phi;
aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aGridSize=GridDensity*OrderOfApproxGamma;
Gamma(1) = fundefn(ApproxMethod,OrderOfApproxGamma ,aMin,aMax);
aGrid=funnode(Gamma(1));
Gamma0=repmat(((aGrid-phi)./(-2*phi)),1,2);
CoeffGamma=ones(OrderOfApproxGamma,sSize);
for inx_s=2:sSize
    Gamma(inx_s) = fundefn(ApproxMethod,OrderOfApproxGamma ,aMin,aMax);
    
    CoeffGamma(:,inx_s)=FitGammaCoeff(Gamma0(:,inx_s),aGrid,Gamma(inx_s),phi,Para);
    CoeffGamma(:,inx_s)=FitGammaCoeffLP(Gamma0(:,inx_s),aGrid,Gamma(inx_s),phi,Para);
end
 
Eqb.CoeffConsumptionPolicy=CoeffConsumptionPolicy;
Eqb.C=C;
Eqb.CoeffAPolicy=CoeffAPolicy;
Eqb.A=A;
Eqb.CoeffGamma=CoeffGamma;
Eqb.Gamma=Gamma;
Eqb.q=q;
Eqb.phi=phi;

[~,Eqb]=ResBondMarketPrice(q,Eqb,Para);
[resINQ,resEQ] =ResBondMarketPriceKnitro(q,Eqb,Para)
Error=ComputeError(Eqb,Para);

DummyObj=@(q) 1;
lb=delta^2;
ub=1;
opts = optimset('Algorithm', 'Interior-point', 'Display','iter','TolCon',1e-6);    
q=delta
[q,fval,exitflag]=ktrlink(@(q) DummyObj(q),q,[],[],[],[],lb,ub,@(q)ResBondMarketPriceKnitro(q,Eqb,Para) ,opts);
options=optimset('Display','iter','TolX',Para.ErrorTol);
%q=fzero(@(q) ResBondMarketPrice(q,Eqb,Para),q,options)
[~,Eqb]=ResBondMarketPrice(q,Eqb,Para);
Error=ComputeError(Eqb,Para);
Error.ConsumptionPolicy
Error.Gamma
Error.APolicy
Error.MarketClearing
CoeffConsumptionPolicy=Eqb.CoeffConsumptionPolicy;
C=Eqb.C;
CoeffAPolicy=Eqb.CoeffAPolicy;
A=Eqb.A;
CoeffGamma=Eqb.CoeffGamma;
Gamma=Eqb.Gamma;
phi=Eqb.phi;
q=Eqb.q;
aMin=phi;
aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aGridSize=GridDensity*OrderOfApproxConsumptionPolicy;
aGrid=linspace(aMin,aMax,aGridSize)';

figure()
for inx_s=1:sSize
    subplot(1,2,inx_s)
plot(aGrid,funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid))
end


figure()
for inx_s=1:sSize
    subplot(1,2,inx_s)
plot(aGrid,funeval(CoeffAPolicy(:,inx_s),A(inx_s),aGrid),'k')
xlabel('a')
ylabel('A(a)')
line([phi;phi],[aMin,aMax],'LineWidth',2)
axis([aMin*2 aMax aMin aMax])
end


figure()
for inx_s=1:sSize
    subplot(1,2,inx_s)
plot(aGrid,funeval(CoeffGamma(:,inx_s),Gamma(inx_s),aGrid))
end
plot(funeval(CoeffGamma(:,inx_s),Gamma(inx_s),aGrid))

% Update q, phi
%% CHECK CONVERGENCE
%% PLOT DELIVERABLES
%% UPDATE
%ErroTotal=1;
%err_tol=1e-3;
%ctr=1;
%adjfactor=0;
%err_tol=1e-1;
% while ErroTotal>err_tol&& ctr <500
% phi=(-S(1)/(1-q))*NonZeroAdj;
% tic
% for ctrIn=1:NumIter
% [C,CoeffConsumptionPolicy]=UpdateConsumptionCoeff(C,CoeffConsumptionPolicy,CoeffAPolicy,A,q,phi,Para);
% [A,CoeffAPolicy]=UpdateACoeff(CoeffConsumptionPolicy,C,A,CoeffAPolicy,q,phi,Para);
% [Gamma,CoeffGamma]=UpdateGammaCoeff(CoeffAPolicy,A,Gamma,CoeffGamma,phi,Para);
% end
% 
% adjfactorold=adjfactor;
% res=ResMarketClearing(CoeffGamma,Gamma, CoeffAPolicy,A,phi,Para);
% adjfactor=(abs(res/phi))*1e-2;
% DirChange=sign(res);
% q=q*grelax+(1-grelax)*(max(min(q*(1+DirChange*adjfactor),1),delta));
% adjfactornew=adjfactor;
% err(ctr)=abs(res);
% ErroTotal=err(ctr);
% ctr=ctr+1;
% toc
% end