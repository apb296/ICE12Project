% This file solves the Huggert model with projection methods 
clc
clear all
close all
warning off all

%global Eqb;
%% SET PARAMETERS 
% Technology and Prefernces
sigma=5; % Risk Aversion
delta=.96; % Time discount factor
S=[1  5]; % Endowment shocks
sSize=length(S);
alpha=.9; % Probability of staying in the same state
P=[alpha 1-alpha;1-alpha alpha]; % Stochastic shock matrix
%P=[.5 .5 ;.5 .5]; % Stochastic shock matrix
%P=[.95 .05 ;.075 1-.075]; % Stochastic shock matrix
Para.sigma=sigma;
Para.delta=delta;
Para.S=S;
Para.sSize=sSize;
Para.P=P;
Para.flagNaturalBorrowingLimit=1;

%Para.AdhocBorrowingLimit=(-S(1)/(1-delta))*.2;

% others 
ApproxMethod='cheb';
OrderOfApproxConsumptionPolicy=6;
%OrderOfApproxAPolicy=7;
OrderOfApproxGamma=6;
ErrorTol=1e-4;
GridDensity=1;
NonZeroAdj=.98;
Para.ApproxMethod=ApproxMethod;
Para.OrderOfApproxConsumptionPolicy=OrderOfApproxConsumptionPolicy;
%Para.OrderOfApproxAPolicy=OrderOfApproxAPolicy;
Para.GridDensity=GridDensity;
Para.OrderOfApproxGamma=OrderOfApproxGamma;
Para.NonZeroAdj=NonZeroAdj;
grelax=.9;
NumIter=100;
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
% 4) Gamma[a,s] - Ergodic distribution over a,s

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
%aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aMax=-phi;
aGridSize=GridDensity*OrderOfApproxConsumptionPolicy;
C(1) = fundefn(ApproxMethod,OrderOfApproxConsumptionPolicy ,aMin,aMax);
aGrid=funnode(C(1));
CoeffConsumptionPolicy=ones(OrderOfApproxConsumptionPolicy,sSize);
C0=repmat((1-q).*aGrid,1,sSize)+repmat(S,aGridSize,1);
CoeffConsumptionPolicy(:,1)=funfitxy(C(1),aGrid,C0(:,1));


for inx_s=2:sSize
    C(inx_s) = fundefn(ApproxMethod,OrderOfApproxConsumptionPolicy ,aMin,aMax);
    CoeffConsumptionPolicy(:,inx_s)=funfitxy(C(inx_s),aGrid,C0(:,inx_s));
end



% -- Gamma[a,s,q]---
aMin=phi;
%aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aMax=-phi;
aGridSize=GridDensity*OrderOfApproxGamma;
Gamma(1) = fundefn(ApproxMethod,OrderOfApproxGamma ,aMin,aMax);
%    aGrid=funnode(Gamma(inx_s));
 %   Gamma0(:,inx_s)=(normcdf(aGrid,0,(aMax-aMin)/2));
    %CoeffGamma(:,inx_s)=FitGammaCoeff(Gamma0(:,inx_s),aGrid,Gamma(inx_s),phi,Para);

aGrid=funnode(Gamma(1));
GammaPDF=repmat(normpdf(aGrid,0,(aMax-aMin)/4),1,2);
GammaPDF=GammaPDF./repmat(sum(GammaPDF),aGridSize,1);

for n=1:aGridSize
Gamma0(n,:)=sum(GammaPDF(1:n,:),1);
end
CoeffGamma=ones(OrderOfApproxGamma,sSize);
CoeffGamma(:,1)=FitGammaCoeffLP(Gamma0(:,1),aGrid,Gamma(1),phi,Para);

%Gamma0=repmat(((aGrid-phi)./(-2*phi)),1,2);

for inx_s=2:sSize
    Gamma(inx_s) = fundefn(ApproxMethod,OrderOfApproxGamma ,aMin,aMax);
%    aGrid=funnode(Gamma(inx_s));
 %   Gamma0(:,inx_s)=(normcdf(aGrid,0,(aMax-aMin)/2));
    %CoeffGamma(:,inx_s)=FitGammaCoeff(Gamma0(:,inx_s),aGrid,Gamma(inx_s),phi,Para);
    CoeffGamma(:,inx_s)=FitGammaCoeffLP(Gamma0(:,inx_s),aGrid,Gamma(inx_s),phi,Para);
end
figure()
plot(funeval(CoeffGamma(:,1),Gamma(1),aGrid),'r')
hold on
plot(funeval(CoeffGamma(:,1),Gamma(1),aGrid,1),'g')
plot(Gamma0(:,1),'k')
hold off

Eqb.CoeffConsumptionPolicy=CoeffConsumptionPolicy;
Eqb.C=C;
Eqb.CoeffGamma=CoeffGamma;
Eqb.Gamma=Gamma;
Eqb.q=q;
Eqb.phi=phi;


% qGrid=linspace(delta,.99,5);
% for i_q=1:5
%     tic
%     resq(i_q)=ResBondMarketPriceBisection(qGrid(i_q),Eqb,Para,'solver');
%     toc
% end
% figure()
% plot(qGrid,resq)
% grid on

%% Iterate

DummyObj=@(q) 1;
lb=delta;
ub=1;
opts = optimset('Algorithm', 'interior-point', 'Display','iter','TolCon',1e-7);    
q=delta;
[q,fval,exitflag]=ktrlink(@(q) DummyObj(q),q,[],[],[],[],lb,ub,@(q)ResBondMarketPriceKnitro(q,Eqb,Para,'solver') ,opts);
%options=optimset('Display','iter','TolX',Para.ErrorTol);
%[q,fval,exitflag]=fzero(@(q)ResBondMarketPriceBisection(q,Eqb,Para,'solver'),[delta],options);
%% SOLUTION
[~,~,Eqb]=ResBondMarketPriceKnitro(q,Eqb,Para,'eqb');
Error=ComputeError(Eqb,Para);
Error.ConsumptionPolicy
Error.Gamma
Error.MarketClearing
CoeffConsumptionPolicy=Eqb.CoeffConsumptionPolicy;
C=Eqb.C;
CoeffGamma=Eqb.CoeffGamma;
Gamma=Eqb.Gamma;
phi=Eqb.phi;
q=Eqb.q;

%% Plots
aMin=phi;
%aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aMax=-phi;
aGridSize=GridDensity*OrderOfApproxConsumptionPolicy;
aGrid=linspace(aMin,aMax,aGridSize)';
figure()
for inx_s=1:sSize
    subplot(1,2,inx_s)
plot(aGrid,funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid),'k','LineWidth',2)
xlabel('a')
ylabel('C[a,s]')
title(['s=' num2str(inx_s)])
end

figure()
for inx_s=1:sSize
    subplot(1,2,inx_s)
ANew(:,inx_s)=(((aGrid+S(inx_s)-funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid))/q)); % Savings given a,s
plot(aGrid,ANew(:,inx_s),'k','LineWidth',2)
hold on
plot(aGrid,aGrid,':k','LineWidth',2)
xlabel('a')
ylabel('A[a,s]')
title(['s=' num2str(inx_s)])
axis tight
end



figure()
for inx_s=1:sSize
    subplot(1,2,inx_s)
ANew(:,inx_s)=(((aGrid+S(inx_s)-funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid))/q)); % Savings given a,s
plot(aGrid,ANew(:,inx_s)-aGrid,'k','LineWidth',2)
xlabel('a')
ylabel('A[a,s]-a')
title(['s=' num2str(inx_s)])
axis tight
end




aMin=phi;
%aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aMax=-phi;
aGridSize=GridDensity*OrderOfApproxGamma;
aGrid=linspace(aMin,aMax,aGridSize)';
GammaPDF=repmat(normpdf(aGrid,0,(aMax-aMin)/2),1,2);
GammaPDF=GammaPDF./repmat(sum(GammaPDF),aGridSize,1);
for n=1:aGridSize
Gamma0(n,:)=sum(GammaPDF(1:n,:),1);
end

figure()
for inx_s=1:sSize
    subplot(1,2,inx_s)
plot(aGrid,funeval(CoeffGamma(:,inx_s),Gamma(inx_s),aGrid))
hold on
plot(aGrid,Gamma0(:,inx_s),':k','LineWidth',2);
xlabel('a')
ylabel('G[a,s]','Interpreter','Latex')
title(['s=' num2str(inx_s)])
axis tight
end

%% SIMULATION
[~,~,Eqb]=ResBondMarketPriceKnitro(q,Eqb,Para,'eqb');
q=Eqb.q;
aMin=Eqb.phi;
%aMax=(Para.S(2)/(1-Para.delta))*Para.NonZeroAdj;
aMax=-Eqb.phi;
NumSim=5000;
AHist=zeros(NumSim,1);
sHist=zeros(NumSim,1);
sHist(1)=1;
AHist(1)=aMin*.5;
sHist_ind=rand(NumSim,1);

for inx_sim=2:NumSim
    s=sHist(inx_sim-1);
    a=AHist(inx_sim-1);
    c=funeval(CoeffConsumptionPolicy(:,s),C(s),a);
    AHist(inx_sim)=max(min((a+S(s)-c)/q,aMax),aMin);
    if sHist_ind(inx_sim)<P(s,1)
        sHist(inx_sim)=1;
    else
        sHist(inx_sim)=2;
    end
end
BurnSample=.2;
[f,x] = ecdf(AHist(NumSim*BurnSample:NumSim));
figure()
plot(AHist)
figure()
plot(S(sHist))
figure()
plot(x,f)


%% COMPARITIVE STATICS
sigmaMin=3;
sigmaMax=7;
sigmaGridSize=10;
sigmaGrid=linspace(sigmaMin,sigmaMax,sigmaGridSize);

for inx_sigma=1:sigmaGridSize
qGuess=Eqb.q;
Para.sigma=sigmaGrid(inx_sigma);
opts = optimset('Algorithm', 'active-set', 'Display','iter','TolCon',1e-5);    
[q,fval,exitflag]=ktrlink(@(q) DummyObj(q),q,[],[],[],[],lb,ub,@(q)ResBondMarketPriceKnitro(q,Eqb,Para,'solver') ,opts);
qSigma(inx_sigma)=q;
[~,~,Eqb]=ResBondMarketPriceKnitro(q,Eqb,Para,'eqb');
end
figure()
plot(sigmaGrid,qSigma,'LineWidth',2)
hold on
plot(sigmaGrid,delta*ones(sigmaGridSize,1),':k','LineWidth',2)
xlabel('$\sigma$','Interpreter','Latex')
ylabel('q','Interpreter','Latex')



% 
% figure()
% for inx_s=1:sSize
%     subplot(1,2,inx_s)
% plot(aGrid,funeval(CoeffGamma(:,inx_s),Gamma(inx_s),aGrid,1))
% hold on
% aMin=phi;
% aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
% aGridSize=GridDensity*OrderOfApproxGamma;
% aGrid=linspace(aMin,aMax,aGridSize)';
% GammaPDF=normpdf(aGrid,0,(aMax-aMin)/3);
% GammaPDF=GammaPDF./sum(GammaPDF);
% plot(aGrid,GammaPDF(:,inx_s),':k','LineWidth',2);
% end
