clc
clear
% options  =  optimoptions( 'fmincon', ...
%                           'FiniteDifferenceStepSize', 1e-7, ...
%                           'OptimalityTolerance', 1e-7, ...
%                           'StepTolerance', 1e-7, ...
%                           'FunctionTolerance', 1e-7, ...
%                           'MaxIter',200, ...
%                           'SpecifyObjectiveGradient', false, ...
%                           'Display', 'on');
 options=optimset('fmincon');
 options.MaxFunEvals =100000000000;
 options.MaxIter =15;
 options.TolFun=1e-10;
% options.StepTolerance=10^-4;
 %options.exitflag=4;
% options.OptimalityTolerance=10^(-5);
 options.Display = 'iter';
 options.Jacobian='off';
%options.Algorithm='levenberg-marquardt';
%options.Algorithm='trust-region-dogleg';
init=19;
[x,f,ex]=fmincon(@(var)Workspace_generation(var),init,[],[],[],[],0,inf,[],options);
x