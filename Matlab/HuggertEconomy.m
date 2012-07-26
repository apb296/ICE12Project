% This file solves the Huggert model with projection methods 
clc
clear all
close all
warning off all

%global Eqb;
%% SET PARAMETERS 
% Technology and Prefernces
sigma=2; % Risk Aversion
delta=.9; % Time discount factor
S=[1  1.5]; % Endowment shocks
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
OrderOfApproxConsumptionPolicy=5;
%OrderOfApproxAPolicy=7;
OrderOfApproxGamma=5;
ErrorTol=1e-6;
GridDensity=50;
NonZeroAdj=.9;
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
aGrid=linspace(aMin,aMax,aGridSize)';
C0=repmat((1-q).*aGrid,1,sSize)+repmat(S,aGridSize,1);
CoeffConsumptionPolicy=ones(OrderOfApproxConsumptionPolicy,sSize);
for inx_s=1:sSize
    C(inx_s) = fundefn(ApproxMethod,OrderOfApproxConsumptionPolicy ,aMin,aMax);
    CoeffConsumptionPolicy(:,inx_s)=funfitxy(C(inx_s),aGrid,C0(:,inx_s));
end



% -- Gamma[a,s,q]---
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
%Gamma0=repmat(((aGrid-phi)./(-2*phi)),1,2);
CoeffGamma=ones(OrderOfApproxGamma,sSize);
for inx_s=1:sSize
    Gamma(inx_s) = fundefn(ApproxMethod,OrderOfApproxGamma ,aMin,aMax);
%    aGrid=funnode(Gamma(inx_s));
 %   Gamma0(:,inx_s)=(normcdf(aGrid,0,(aMax-aMin)/2));
    CoeffGamma(:,inx_s)=FitGammaCoeff(Gamma0(:,inx_s),aGrid,Gamma(inx_s),phi,Para);
    %CoeffGamma(:,inx_s)=FitGammaCoeffLP(Gamma0(:,inx_s),aGrid,Gamma(inx_s),phi,Para);
end
figure()
plot(funeval(CoeffGamma(:,1),Gamma(1),aGrid),'r')
hold on
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
opts = optimset('Algorithm', 'active-set', 'Display','iter','TolCon',1e-5);    
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


%% Euler- Lagrange -Philipp  Equation Errors

aMin=Eqb.phi;
%aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aMax=-Eqb.phi;
aGridSize=GridDensity*OrderOfApproxConsumptionPolicy;
aGrid=aMin+(aMax-aMin)*rand(aGridSize,1);
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
   EulerEquationErrors(:,inx_s)=CNew(:,inx_s)-funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid);
end

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
xlabel('a','FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
ylabel('C[a,s]','FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
title(['s=' num2str(inx_s)],'FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
axis tight
end
print(gcf,'-dpng',[ plotpath 'FigConsumptionPolicy.png'])
figure()
for inx_s=1:sSize
    subplot(1,2,inx_s)
ANew(:,inx_s)=(((aGrid+S(inx_s)-funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid))/q)); % Savings given a,s
plot(aGrid,ANew(:,inx_s),'k','LineWidth',2)
hold on
plot(aGrid,aGrid,':k','LineWidth',2)
xlabel('a','FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
ylabel('A[a,s]','FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
title(['s=' num2str(inx_s)],'FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
axis tight
end

print(gcf,'-dpng',[ plotpath 'FigSavingsPolicy.png'])


figure()
for inx_s=1:sSize
    subplot(1,2,inx_s)
ANew(:,inx_s)=(((aGrid+S(inx_s)-funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid))/q)); % Savings given a,s
plot(aGrid,ANew(:,inx_s)-aGrid,'k','LineWidth',2)
xlabel('a','FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
ylabel('A[a,s]-a','FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
title(['s=' num2str(inx_s)],'FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
axis tight
end
print(gcf,'-dpng',[ plotpath 'FigDeltaSavings.png'])



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
xlabel('a','FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
ylabel('G[a,s]','Interpreter','Latex','FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
title(['s=' num2str(inx_s)],'FontWeight','bold','FontSize',12,...
    'FontName','Bitstream Charter')
axis tight
end
print(gcf,'-dpng',[ plotpath 'FigGamma.png'])
%% SIMULATION
[~,~,Eqb]=ResBondMarketPriceKnitro(q,Eqb,Para,'eqb');
q=Eqb.q;
aMin=Eqb.phi;
%aMax=(Para.S(2)/(1-Para.delta))*Para.NonZeroAdj;
aMax=-Eqb.phi;
NumSim=15000;
AHist=zeros(NumSim,1);
sHist=zeros(NumSim,1);
sHist(1)=1;
AHist(1)=0;
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
BurnSample=.5;
[f,x] = ecdf(AHist(NumSim*BurnSample:NumSim));
% figure()
T=100;
%plot(AHist(end-T:end));
X.Data=AHist(end-T:end);
    X.sHist=sHist(end-T:end);
    X.name={'$a$'};
    PlotSimul(X)
    title('Assets');
print(gcf,'-dpng',[ plotpath 'FigSimulationA.png'])
% figure()
% plot(S(sHist))
figure()
plot(x,f,':k','LineWidth',2)
hold on
plot(aGrid,funeval(CoeffGamma(:,1),Gamma(1),aGrid)*.5+funeval(CoeffGamma(:,2),Gamma(2),aGrid)*.5,'k','LineWidth',2)
axis tight
print(gcf,'-dpng',[ plotpath 'FigSimulationCDF.png'])

%% COMPARITIVE STATICS
sigmaMin=sigma*.5;
sigmaMax=sigma*2;
sigmaGridSize=5;
sigmaGrid=linspace(sigmaMin,sigmaMax,sigmaGridSize);

for inx_sigma=1:sigmaGridSize
qGuess=Eqb.q;
Para.sigma=sigmaGrid(inx_sigma);
opts = optimset('Algorithm', 'active-set', 'Display','iter','TolCon',1e-5);    
[q,fval,exitflag]=ktrlink(@(q) DummyObj(q),q,[],[],[],[],lb,ub,@(q)ResBondMarketPriceKnitro(q,Eqb,Para,'solver') ,opts);
qSigma(inx_sigma)=q;
[~,~,Eqb]=ResBondMarketPriceKnitro(q,Eqb,Para,'eqb');
[Error]=ComputeError(Eqb,Para);
L2ErrorConsumptionPolicy(inx_sigma)=Error.ConsumptionPolicy;
L2ErrorGammaPolicy(inx_sigma)=Error.Gamma;
L2MarketClearingError(inx_sigma)=Error.MarketClearing;
CoeffConsumptionPolicy=Eqb.CoeffConsumptionPolicy;
C=Eqb.C;
CoeffGamma=Eqb.CoeffGamma;
Gamma=Eqb.Gamma;
phi=Eqb.phi;
q=Eqb.q;


aMin=Eqb.phi;
%aMax=(S(2)/(1-delta))*Para.NonZeroAdj;
aMax=-Eqb.phi;
aGridSize=GridDensity*OrderOfApproxConsumptionPolicy;
aGrid=aMin+(aMax-aMin)*rand(aGridSize,1);
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
   EulerEquationErrors(:,inx_s)=CNew(:,inx_s)-funeval(CoeffConsumptionPolicy(:,inx_s),C(inx_s),aGrid);
end



L2EulerEquationErrors(inx_sigma)=(sum(sqrt(sum(EulerEquationErrors.^2)))/aGridSize);
end


   matrix = [1.5 1.764; 3.523 0.2];
   rowLabels = {'row 1', 'row 2'};
   columnLabels = {'col 1', 'col 2'};
   matrix2latex(matrix, 'out.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');

   
   
ErrorMatrix=[L2ErrorConsumptionPolicy;L2ErrorGammaPolicy;L2MarketClearingError;L2EulerEquationErrors] ;
rowLabels = {'$\|\mathcal{C}^{j+1}-\mathcal{C}^{j}\|_2$','$\|\Gamma^{j+1}-\Gamma^{j}\|_2$','$\frac{1}{\phi} \|\int{d\Gamma[a,s]\mathcal{A}[a,s]\|}$', 'Euler Eq Error'};
columnLabels = {['$\sigma=$' num2str(sigmaGrid(1))] , ['$\sigma=$' num2str(sigmaGrid(2))],['$\sigma=$' num2str(sigmaGrid(3))],['$\sigma=$' num2str(sigmaGrid(4))],['$\sigma=$' num2str(sigmaGrid(5))]};
%matrix2latex([L2ErrorConsumptionPolicy;L2ErrorGammaPolicy;L2MarketClearingError;L2EulerEquationErrors], [texpath 'Diagnostics.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
matrix2latex(ErrorMatrix, [texpath 'Diagnostics.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%6.2e', 'size', 'footnotesize');




figure()
plot(sigmaGrid,qSigma,'LineWidth',2)
%hold on
%plot(sigmaGrid,delta*ones(sigmaGridSize,1),':k','LineWidth',2)
xlabel('$\sigma$','Interpreter','Latex')
ylabel('q','Interpreter','Latex')

print(gcf,'-dpng',[ plotpath 'FigComparitiveStatics.png'])

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
ParameterMatrix=[sigma ; delta ; alpha ];
rowLabels = {'$\sigma$', '$\delta$', '$\alpha$'};
columnLabels = {'Values'};
matrix2latex(ParameterMatrix, [texpath 'Parameters.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'normal');

   
