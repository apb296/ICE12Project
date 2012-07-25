function [CoeffConsumptionPolicy,C]=UpdateConsumptionCoeffLP(C,CoeffConsumptionPolicy,CoeffAPolicy,A,q,phi,Para);
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


Degree=C(1).n; % n

% Choose the points to fit the value and the shape constraints
% Take Chebychef nodes for both of them


FitPoints = aGrid(1:2:Degree);
ShapeTestPoints = aGrid(2:2:Degree);
NumFitPoints=length(FitPoints); % m
NumShapePoints=length(ShapeTestPoints);



%FitPoints=aMin+rand(NumFitPoints,1)*(aMax-aMin);
%ShapeTestPoints=aMin+rand(NumShapePoints,1)*(aMax-aMin);

% Define the linear obj
Penalty = zeros(Degree,1);
Penalty(NumFitPoints+1:Degree)=((NumFitPoints+1:Degree)+1-NumFitPoints).^2';
Obj=ones(Degree,1)+Penalty;

%Equality constraints the first are all set to Gamma0 the second to [0 1]
Aeq = funbas(C(1),FitPoints,0);


beq1 = funeval(CoeffConsumptionPolicy(:,1),C(1),aGrid(1:2:Degree),0);
beq2 = funeval(CoeffConsumptionPolicy(:,2),C(1),aGrid(1:2:Degree),0);

%Shape constraints Aineq <= 0
epsilon=1e-6;
Aineq =-funbas(C(1),ShapeTestPoints,1);
bineq = -epsilon*ones(NumShapePoints,1);


options=optimset('Display', 'off');
CoeffConsumptionPolicy(:,1) = linprog(Obj,Aineq,bineq,Aeq,beq1,[],[],[],options);
CoeffConsumptionPolicy(:,2) = linprog(Obj,Aineq,bineq,Aeq,beq2,[],[],[],options);




end

