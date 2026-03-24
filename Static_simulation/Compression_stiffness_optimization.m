function r=Workspace_generation(S)
global sig1 
global stab

Param.E = 3*10^9; % Young modulu
Param.G=8.543*10^4; % shear modulu
%Param.r = 0.02; % rod radius
Param.rho = 1240*0.2; % Mass density

%Platform parameters
Param.a=0.02;
Param.b=0.1;
Param.c=0.02;
Param.rhop=1.41*10^3;
Param.Vp=Param.a*Param.b*Param.c;
Param.mp=Param.rhop*Param.Vp;
Param.Jp1=1/12*Param.mp*(Param.b^2+Param.c^2);
Param.Jp2=1/12*Param.mp*(Param.a^2+Param.c^2);
Param.Jp3=1/12*Param.mp*(Param.a^2+Param.b^2);
Param.Mp =diag([Param.Jp1,Param.Jp2,Param.Jp3,Param.mp,Param.mp,Param.mp]);%cross-sectional inertia


Param.g =[0 ; 0 ; 0 ; 0 ; 0; -9.81]; % gravity
Param.L = 0.2; % Rod length
Param.d=0.001;
Param.A = Param.c*Param.d; % Cross section area
Param.J1 = Param.c*Param.d*(Param.c^2+Param.d^2)/12; % Polar inertia moment
Param.J2 = Param.c^3*Param.d/12; % Inertia moment
Param.J3 = Param.c*Param.d^3/12; % Inertia moment
%Param.mu=1000;
%Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 ,S, Param.G*Param.A , Param.G*Param.A]);
Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia
Param.dX=Param.L/100; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);



%Param.Fep=[0;Param.mp*9.81;0];
%Param.Fep=0*[0;40;0];
%Param.mu=10^-3;
%Tendons parameters
Param.Rb= 0.025;
Param.c1=0*-pi/8;
Param.c2=0*pi/8;
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)]; % coordinate of the tendon in a cross-section
Param.n=3;  %strain modes number
Param.na=2; %number of actuated strains
Param.dX=Param.L/100;%length of segments
Param.n_lambda=6; %dimension of Lambda
Y=0:Param.dX:Param.L;
Param.Tendons_list=zeros(3*(Param.n_X+1),2);
for i=1:Param.n_X+1
Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi,Param.Rb)]; % Parallel case
%Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-i*Param.dX))]; % Convergent routing
%Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-2*i*Param.dX))]; % Convergent routing
%Parallel truncated
%if (i-1)*Param.dX<Param.L/2
%    Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi/2,Param.Rb),Param.Tendon_coordinate(pi,Param.Rb),Param.Tendon_coordinate(3*pi/2,Param.Rb)];
%end
end
  
%Param.Forces_Tendons=[0
%5
%5
%0]; 
   
Param.Ftip=0*[0;-1/sqrt(2);1/sqrt(2);0;0;0];
%Param.Forces_Tendons(1)=1;


%Initial conditions
Param.theta01=0*-pi/8;
Param.theta02=0*pi/8;
R01=[cos(Param.theta01) -sin(Param.theta01) 0; sin(Param.theta01) cos(Param.theta01) 0; 0 0 1];
R02=[cos(Param.theta02) -sin(Param.theta02) 0; sin(Param.theta02) cos(Param.theta02) 0; 0 0 1];
Param.g1=eye(4);
Param.g1(1:3,1:3)=R01;
Param.g2=eye(4);
Param.y2=0.1;
Param.g2(1:3,1:3)=R02;
Param.g2(2,4)=Param.y2;
Param.CM1=[Param.a/2;Param.y2/2];
Param.CM2=[Param.a/2;-Param.y2/2];
%Param.B =[eye(3);zeros(3)]; 
%Param.B_bar = [zeros(3);eye(3)];
% Param.B     = [0 0 0;1 0 0;0 0 0;0 1 0;0 0 0;0 0 1];
% Param.B_bar = [1 0 0;0 0 0;0 1 0;0 0 0;0 0 1;0 0 0];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;0];

% Param.B     = [0;1;0;0;0;0];
% Param.B_bar = [1 0 0 0 0;0 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;1;0;0];

Param.B     = [0 0;0 0;1 0;0 1;0 0;0 0];
Param.B_bar = [1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0;0 0 1 0;0 0 0 1];
Param.xi_0=[0;0;0;1;0;0];
Param.xi_a0=Param.B'*Param.xi_0;
Param.xi_c=[0;0;0;0];



Param.Ha=Param.B'*Param.H*Param.B;%the matrix of the reduced Hooke coefficients
%Param.Da=Param.B'*Param.D*Param.B;%reduced damping matrix
Y=0:Param.dX:Param.L;
SK=zeros(Param.n*Param.na,Param.n*Param.na);
for k=1:Param.n_X+1
    phival=Phi(Param.na,Param.n,Y(k),Param.L);
    fK=phival'*Param.Ha*phival;
    if k==1
      %  fK10=fK1;
        fK0=fK;
    elseif k==(Param.n_X+1)
       % fK1n=fK1;
        fKn=fK;
    else
       % SK1=SK1+fK1;
        SK=SK+fK;
    end
end
Param.Keps= Param.dX.*((fK0+fKn)./2 + SK);
%Trapeze integration method
toc

lb=zeros(2*Param.n*Param.na+Param.n_lambda+4,1);
ub=zeros(2*Param.n*Param.na+Param.n_lambda+4,1);
lb(1:2*Param.n*Param.na+Param.n_lambda)=-inf;
lb(2*Param.n*Param.na+Param.n_lambda+1:2*Param.n*Param.na+Param.n_lambda+4)=0;
ub(1:2*Param.n*Param.na+Param.n_lambda)=inf;
ub(2*Param.n*Param.na+Param.n_lambda+1:2*Param.n*Param.na+Param.n_lambda+4)=100;
%init=zeros(19,1);
Param.mres=0;
Param.T_min=[0;0];
Param.T_max=[100;100];
%init(:,1)=[zeros(2*Param.n*Param.na+Param.n_lambda+3,1);sqrt(Param.T_min);sqrt(Param.T_max)];
%init(:,1)=[zeros(2*Param.n*Param.na+Param.n_lambda+1,1);Param.T_min;Param.T_max];
%init(:,1)=[zeros(2*Param.n*Param.na+Param.n_lambda+3,1);0.00001;0.00001;Param.T_max];
%init(:,1)=[zeros(2*Param.n*Param.na+Param.n_lambda+3,1);0.5;0.5;10;10];
%init(:,1)=zeros(2*Param.n*Param.na+Param.n_lambda+3,1);%;0;0;sqrt(Param.T_max);sqrt(Param.T_max)]];
% init=[ -1.504989513688684 %T=4N
  % -0.000000000000000
  % 18.155268566737533
  % -0.000132577106367
  % -0.000000000000000
  %  0.000001406257921
  % -1.504989513688684
  % -0.000000000000001
  % -2.332743627614550
  %  0.000069442963581
  %  0.000000000000000
  % -0.000000409061144
  % -0.300997902737737
  %  0.201877600262403
  % -0.015130306863276
  %  0.081143386943996
  % -4.127537878744650
  %  0.625922978500879
  % -0.010975243005337
  %  4.127537878744650
  % -0.625922978500879];

 %     init=[ -10.104343199150291 %T=4N
 %   0.000000000463339
 %  27.457382273523919
 %  -0.000095091781586
 %  -0.000000000000000
 %   0.000019186549886
 % -10.104343199150291
 %   0.000000000337865
 % -11.868133473848227
 %   0.000041766350546
 %   0.000000000000000
 %  -0.000012372881256
 %  -2.020868639830058
 %   4
 %   0
 %   0.083392714504166
 %  -1.478386599850913
 %   2.356167496296465
 %  -0.058260142305702
 %   1.478386599838229
 %  -2.356167496297970];
init(:,1)=zeros(2*Param.n*Param.na+Param.n_lambda+3,1);%0.5;0.5;10;10];
init(2*Param.n*Param.na+2,1)=0.2+Param.a/2;
init(2*Param.n*Param.na+3,1)=Param.y2/2;
% P_des=[0.2+Param.a/2;Param.y2/2];
% plot(0.2+Param.a/2,Param.y2/2,'r.'); title('Workspace');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
% xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;


figure(1)
hold on
plot(0.2+Param.a/2,Param.y2/2,'b.');   
%plot(                0.076117184007060,-0.120248612677159 ,'b.');
title('Type2');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;

figure(2)
hold on
plot(0.2+Param.a/2,Param.y2/2,'b.');   
%plot(            0.076117184007060,-0.120248612677159,'b.');
title('Stability');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;

%Ws(:,1)=[0.076117184007060,-0.120248612677159];
Ws(:,1)=[0.2+Param.a/2,Param.y2/2];
s=0.01;
k=1;
i=1;
while i<=length(Ws(1,:))
    Neighbors=[Ws(:,i)+[-s;0],Ws(:,i)+[s;0],Ws(:,i)+[0;s],Ws(:,i)+[0;-s],Ws(:,i)+[-s;s],Ws(:,i)+[-s;-s],Ws(:,i)+[s;-s],Ws(:,i)+[s;s]];%,Ws(:,i+1)+[0;-0.005]];
    %k=k-1;
    j=1;
    while j<9
    %Je fais tous les voisins et je mets ceux qui sont dans le workspace
    %dans le vecteur WS. Je jette le reste ensuite prends comme CI le
    %prochain point du vecteur WS et ainsi de suite.
    t=1;
    for l=1:length(Ws(1,:))
        if abs(Neighbors(:,j)-Ws(:,l))<[10^(-10);10^(-10)]
            t=0;
        end
    end
 if t==1 & Neighbors(2,j)<=Param.y2/2
 tic
 P_des=Neighbors(:,j);
 options=optimset('fsolve');
 options.MaxFunEvals =100000000000;
 options.MaxIter =15;
 options.TolFun=1e-10;
% options.StepTolerance=10^-4;
 %options.exitflag=4;
% options.OptimalityTolerance=10^(-5);
 options.Display = 'iter';
 options.Jacobian='on';
%options.Algorithm='levenberg-marquardt';
%options.Algorithm='trust-region-dogleg';
%sing=1;
% options  =  optimoptions( 'fsolve', ...
%                           'FiniteDifferenceStepSize', 1e-7, ...
%                           'OptimalityTolerance', 1e-7, ...
%                           'StepTolerance', 1e-7, ...
%                           'FunctionTolerance', 1e-7, ...
%                           'MaxIter',200, ...
%                           'SpecifyObjectiveGradient', true, ...
%                           'Display', 'off');
[x,f,a,b,J1]=fsolve(@(var)Static_Gradient_Limit_actuation_V3(var,P_des,Param),init(:,i),options);
%x;
%[x,f,a,b,J1]=fsolve(@(var)Static_Gradient_Limit_actuation(var,P_des,Param),init(:,i),options);
%[x,f,a,b,J1]=fsolve(@(var)Static_Gradient(var,P_des,Param),init(:,i),options);
%options.Jacobian='off';
%[x,f,a,b,J2]=fsolve(@(var)Static(var,P_des,Param),init(:,i),options);
toc
%J=J1-J2
%[x,f,a,b,J2]=fsolve(@(var)Static_Limited_actuation(var,P_des,Param),init(:,i),options);
%[U1,SS1,V1]=svd(J1);
%[U2,SS2,V2]=svd(J2);
%D1=det(J1)
%D2=det(J2)
%e1=U1-U2;
%e2=SS1-SS2;
%e3=V1-V2;
%gam1=eig(J1)
%gam2=eig(J2)
%max(abs(f))
%invj1=inv(J1)
%invJ=J2*inv(J1)
x
q1=x(1:Param.n*Param.na);
q2=x(Param.n*Param.na+1:2*Param.n*Param.na);
q=[q1;q2];
qp=x(2*Param.n*Param.na+1);
T=x(2*Param.n*Param.na+2:2*Param.n*Param.na+3);
Lambda=x(2*Param.n*Param.na+4:2*Param.n*Param.na+3+Param.n_lambda);

% S1=x(2*Param.n*Param.na+1+Param.n_lambda+1:2*Param.n*Param.na+1+Param.n_lambda+2);
% S2=x(2*Param.n*Param.na+1+Param.n_lambda+3:2*Param.n*Param.na+1+Param.n_lambda+4);
    if max(abs(f))<1e-4 & T(1)>=0 & T(2)>=0 & T(1)<=50 & T(2)<=50
       % hold on
        %             rx=(r1(1,end)+r2(1,end))/2;
        %             ry=(r1(2,end)+r2(2,end))/2;
%         if sing==1
%         plot(P_des(1),P_des(2),'r.'); title('Workspace');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
%         xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;
%         elseif sing ==2
%           plot(P_des(1),P_des(2),'g.'); title('Workspace');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
%         xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;
%        
%         else
%          plot(P_des(1),P_des(2),'b.'); title('Workspace');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
%          xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;   
%          end
        %c=[10^(-10):10^(-10):10^(-5)];
% sig1
% figure(1)
%         scatter(T(1),T(2),[],sig1,'filled')
%         %        scatter(P_des(1),P_des(2),[],T(1),'filled')
%         colorbar
%         colormap jet;
%         drawnow;
% figure(1)
%         scatter(P_des(1),P_des(2),[],sig1,'filled')
%         colorbar
%         colormap jet;
%          %caxis([0 10^(-4)]);
%         %caxis([-10^(-5) 10^(-5)])
%         drawnow;
% 
%  figure(2)
%         scatter(P_des(1),P_des(2),[],stab,'filled')
%         %        scatter(P_des(1),P_des(2),[],T(1),'filled')
%         colorbar
%         colormap jet;
%         %caxis([-10^(-5) 10^(-5)])
%         drawnow;

% figure (2)
% %hold on
% % if stab>10^-6
% %     plot(P_des(1),P_des(2),'b.','LineWidth', 4);title('Stability');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
% %          xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;
% % elseif stab<-10^(-6)
% %     plot(P_des(1),P_des(2),'g.','LineWidth', 4);title('Stability');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
% %          xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;
% % else
% %   plot(P_des(1),P_des(2),'r.','LineWidth', 4);title('Stability');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
% %          xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;  
% % end
%         scatter(P_des(1),P_des(2),[],stab,'filled')
%         colorbar
%         colormap jet;
%         caxis([-10^(-10) 10^(-10)])
%         drawnow;
        
% r1=zeros(3,Param.n_X+1);
% r2=zeros(3,Param.n_X+1);
% r2(2,1)=Param.y2;
% Q1=zeros(4,Param.n_X+1);
% Q2=zeros(4,Param.n_X+1);
% Q1(:,1)=[cos(Param.theta01/2);0;0;sin(Param.theta01/2)];
% Q2(:,1)=[cos(Param.theta02/2);0;0;sin(Param.theta02/2)];
% for p=1:Param.n_X
%   %   g=Geometric_model((j-1)*Param.dX,q,Param);
%   %   r(:,j)=g(1:3,4);
%         phi=Phi(Param.na,Param.n,(p-1)*Param.dX,Param.L);%Functions basis values at X
%         xia1=Param.xi_a0 + phi*q1;
%         xi1=Param.B*xia1+Param.B_bar*Param.xi_c;
%         xia2=Param.xi_a0 + phi*q2;
%         xi2=Param.B*xia2+Param.B_bar*Param.xi_c;
%         %xi(4:6)=[1;0;0];
%         K1=xi1(1:3); %angular strain
%         Gamma1=xi1(4:6);%Linear strain 
%         K2=xi2(1:3); %angular strain
%         Gamma2=xi2(4:6);%Linear strain 
% 
%         R1=eye(3) + 2/(Q1(:,p)'*Q1(:,p)) * [-Q1(3,p)^2-Q1(4,p)^2, Q1(2,p)*Q1(3,p)-Q1(4,p)*Q1(1,p),Q1(2,p)*Q1(4,p) + Q1(3,p)*Q1(1,p) ; 
%                                  Q1(2,p)*Q1(3,p)+Q1(4,p)*Q1(1,p), -Q1(2,p)^2-Q1(4,p)^2,Q1(3,p)*Q1(4,p) - Q1(2,p)*Q1(1,p) ; 
%                                  Q1(2,p)*Q1(4,p)-Q1(3,p)*Q1(1,p), Q1(3,p)*Q1(4,p) + Q1(2,p)*Q1(1,p), -Q1(2,p)^2-Q1(3,p)^2];
%         Q_X1 = [ 0, -K1(1), -K1(2), -K1(3);
%                  K1(1), 0, K1(3), -K1(2);
%                 K1(2), -K1(3), 0, K1(1);
%                 K1(3), K1(2), -K1(1), 0 ] * Q1(:,p)/2;
%         r_X1 = R1*Gamma1;
%         Q1(:,p+1)=Q1(:,p) + Q_X1*Param.dX;
%         r1(:,p+1)=r1(:,p) + r_X1*Param.dX;
%         R1=eye(3) + 2/(Q1(:,p+1)'*Q1(:,p+1)) * [-Q1(3,p+1)^2-Q1(4,p+1)^2, Q1(2,p+1)*Q1(3,p+1)-Q1(4,p+1)*Q1(1,p+1),Q1(2,p+1)*Q1(4,p+1) + Q1(3,p+1)*Q1(1,p+1) ; 
%                                  Q1(2,p+1)*Q1(3,p+1)+Q1(4,p+1)*Q1(1,p+1), -Q1(2,p+1)^2-Q1(4,p+1)^2,Q1(3,p+1)*Q1(4,p+1) - Q1(2,p+1)*Q1(1,p+1) ; 
%                                  Q1(2,p+1)*Q1(4,p+1)-Q1(3,p+1)*Q1(1,p+1), Q1(3,p+1)*Q1(4,p+1) + Q1(2,p+1)*Q1(1,p+1), -Q1(2,p+1)^2-Q1(3,p+1)^2];
% 
% 
%         R2=eye(3) + 2/(Q2(:,p)'*Q2(:,p)) * [-Q2(3,p)^2-Q2(4,p)^2, Q2(2,p)*Q2(3,p)-Q2(4,p)*Q2(1,p),Q2(2,p)*Q2(4,p) + Q2(3,p)*Q2(1,p) ; 
%                                  Q2(2,p)*Q2(3,p)+Q2(4,p)*Q2(1,p), -Q2(2,p)^2-Q2(4,p)^2,Q2(3,p)*Q2(4,p) - Q2(2,p)*Q2(1,p) ; 
%                                  Q2(2,p)*Q2(4,p)-Q2(3,p)*Q2(1,p), Q2(3,p)*Q2(4,p) + Q2(2,p)*Q2(1,p), -Q2(2,p)^2-Q2(3,p)^2];
%         Q_X2 = [ 0, -K2(1), -K2(2), -K2(3);
%                  K2(1), 0, K2(3), -K2(2);
%                 K2(2), -K2(3), 0, K2(1);
%                 K2(3), K2(2), -K2(1), 0 ] * Q2(:,p)/2;
%         r_X2 = R2*Gamma2;
%         Q2(:,p+1)=Q2(:,p) + Q_X2*Param.dX;
%         r2(:,p+1)=r2(:,p) + r_X2*Param.dX;
%         R2=eye(3) + 2/(Q2(:,p+1)'*Q2(:,p+1)) * [-Q2(3,p+1)^2-Q2(4,p+1)^2, Q2(2,p+1)*Q2(3,p+1)-Q2(4,p+1)*Q2(1,p+1),Q2(2,p+1)*Q2(4,p+1) + Q2(3,p+1)*Q2(1,p+1) ; 
%                                  Q2(2,p+1)*Q2(3,p+1)+Q2(4,p+1)*Q2(1,p+1), -Q2(2,p+1)^2-Q2(4,p+1)^2,Q2(3,p+1)*Q2(4,p+1) - Q2(2,p+1)*Q2(1,p+1) ; 
%                                  Q2(2,p+1)*Q2(4,p+1)-Q2(3,p+1)*Q2(1,p+1), Q2(3,p+1)*Q2(4,p+1) + Q2(2,p+1)*Q2(1,p+1), -Q2(2,p+1)^2-Q2(3,p+1)^2];
% end  
% g1=[R1 r1(:,end);0 0 0 1];
% g2=[R2 r2(:,end);0 0 0 1];
% Rp=[cos(qp(1)) -sin(qp(1)); sin(qp(1)) cos(qp(1))];
%     l_rec=0.02;
%     qp=[qp(1);P_des];
%     figure(3)
% 
%     plot(r1(1,:),r1(2,:),'b','LineWidth', 4); title('Planar TACR');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
%     xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
%     hold on
%     plot(r2(1,:),r2(2,:),'b','LineWidth', 4);
%         plot(P_des(1),P_des(2),'g*','LineWidth', 4);
% 
% %       xrec=[qp(2); qp(2)+2*Param.a/2*cos(qp(1))];
% %       yrec=[qp(3); qp(3)+2*Param.a/2*sin(qp(1))];
% %       plot(xrec,yrec,'r','LineWidth', 4);
%      xrec=[qp(2)+(Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1))); qp(2)+(Param.a/2*cos(qp(1))+Param.b/2*sin(qp(1))) ;qp(2)-(Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1)));qp(2)+(-Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1)));qp(2)+(Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1)))];
%      yrec=[qp(3)+(+Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1))) ;qp(3)+(Param.a/2*sin(qp(1))-Param.b/2*cos(qp(1))) ;qp(3)-(Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1)));qp(3)+(-Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1)));qp(3)+(Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1)))];
%      plot(xrec,yrec,'r','LineWidth', 4);
%      hold off

         Ws(:,k+1)=P_des;
         init(:,k+1)=x;
         k=k+1;
        %Neighbors=[Ws(:,i)+[0.005;0],Ws(:,i)+[-0.005;0],Ws(:,i)+[0;0.005],Ws(:,i)+[0;-0.005],Ws(:,i)+[0.005;0.005],Ws(:,i)+[-0.005;0.005],Ws(:,i)+[-0.005;0.005],Ws(:,i)+[-0.005;-0.005]];
        
        %Choisir le plus proche au point trouvé
        %Param.P_des=Neighbors()
        
    end
    end
    j=j+1;
   % pause
    end
%     Param.P_des=Ws(:,i+1);
i=i+1;
end
r=length(Ws(1,:));
end