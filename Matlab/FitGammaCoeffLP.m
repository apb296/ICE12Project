function CoeffGamma=FitGammaCoeffLP(Gamma0,aGrid,Gamma,phi,Para)
aMin=phi;
aMax=(Para.S(2)/(1-Para.delta))*Para.NonZeroAdj;

Degree=Gamma.n; % n

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
Aeq1 = funbas(Gamma,FitPoints,0);
Aeq2 = funbas(Gamma,[aMin; aMax],0);
Aeq = [Aeq1 ; Aeq2];
beq = [Gamma0(1:NumFitPoints);0;1];

%Shape constraints Aineq <= 0
Aineq =-funbas(Gamma,ShapeTestPoints,1);
bineq = zeros(length(ShapeTestPoints),1);

CoeffGamma = linprog(Obj,Aineq,bineq,Aeq,beq);
end
