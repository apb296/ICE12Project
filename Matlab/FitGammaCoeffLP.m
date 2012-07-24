function CoeffGamma=FitGammaCoeffLP(Gamma0,aGrid,Gamma,phi,Para)
aMin=phi;
aMax=(Para.S(2)/(1-Para.delta))*Para.NonZeroAdj;
[CoeffGamma0]=funfitxy(Gamma,aGrid,Gamma0);
Degree=Gamma.n; % n
NumFitPoints=round(.5*Gamma.n)+2; % m
NumShapePoints=Degree;
% Choose the points tot fit the value and the shape constraints
FitPoints=aMin+rand(NumFitPoints,1)*(aMax-aMin);
ShapeTestPoints=aMin+rand(NumShapePoints,1)*(aMax-aMin);

% Define the linear obj
Penalty=((NumFitPoints+1:Degree)+1-NumFitPoints).^2';
Obj=@(CoeffGamma) sum(CoeffGamma(1:NumFitPoints))+...
    sum(Penalty.*(CoeffGamma(NumFitPoints+1:end)));



ShapeConstraints=funbas(Gamma,ShapeTestPoints,1);
LowerBoundryConstraint=funbas(Gamma,aMin);
UpperBoundaryConstraint=funbas(Gamma,aMax);
FitConstraints=funbas(Gamma,FitPoints);
A=-ShapeConstraints;
b=zeros(NumShapePoints,1);
Aeq(1,:)=UpperBoundaryConstraint;
Aeq(2,:)=LowerBoundryConstraint;
Aeq(3:3+NumFitPoints-1,:)=FitConstraints;
beq=[1 ;0;FitPoints];
opts = optimset('Display','off','TolCon',1e-6);    
CoeffGamma = ktrlink(@(CoeffGamma) Obj(CoeffGamma),CoeffGamma0,A,b,Aeq,beq,[],[],[],opts);
end
