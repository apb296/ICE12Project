function CoeffGamma=FitGammaCoeff(Gamma0,aGrid,Gamma,phi,Para)
aMin=phi;
aMax=(Para.S(2)/(1-Para.delta))*Para.NonZeroAdj;
[CoeffGamma0]=funfitxy(Gamma,aGrid,Gamma0);
ShapeTestPoints=aGrid(1:2:end);
%TotalPoints=length(FitPoints)+length(ShapeTestPoints);
%ObjShape=@(CoeffGamma)
%%sum(CoeffGamma(1:FitPoints))+(Degree:NumPoints).^2.*CoeffGamma(Degree:NumPoints);s
LSRes=@(CoeffGamma) (-sum((funeval(CoeffGamma,Gamma,aGrid)-Gamma0).^2))^.5;
n=Gamma.n;
% Test points for shape constraints
%ShapeTestPoints=aMin+rand(NumShapeTestPoints,1)*(aMax-aMin);

NumShapeTestPoints=length(ShapeTestPoints);
ShapeConstraints=funbas(Gamma,ShapeTestPoints,1);
LowerBoundryConstraint=funbas(Gamma,aMin);
UpperBoundaryConstraint=funbas(Gamma,aMax);
A=-ShapeConstraints;
b=zeros(NumShapeTestPoints,1);
%A=[];
%b=[];
Aeq(1,:)=UpperBoundaryConstraint;
Aeq(2,:)=LowerBoundryConstraint;
beq=[1 ;0];
opts = optimset('Algorithm', 'interior-point', 'Display','off','TolCon',1e-6);    
CoeffGamma = ktrlink(@(CoeffGamma) LSRes(CoeffGamma),CoeffGamma0,A,b,Aeq,beq,[],[],[],opts);
end
