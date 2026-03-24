%Main
%Parameters
clc
clear
tic
format long
global sig1 
global sig2
global sig3
%global sing
 %Param.E = 10^10; % Young modulu
% Param.r = 0.005; % rod radius
% Param.g =0*[0 ; 0 ; 0 ; 0 ; 0; -9.81]; % gravity
% Param.L = 1; % Rod length
% Param.A =  pi*Param.r^2; % Cross section area
% Param.J1 = pi*Param.r^4/2; % Polar inertia moment
% Param.J2 = pi*Param.r^4/4; % Inertia moment
% Param.J3 = pi*Param.r^4/4; % Inertia moment
% %Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
% Param.H = diag([1.5708e-2 , 0.7854e-2 , 0.7854e-2 , 3.14e2 ,3.14e2 , 3.14e2]); %Hooke tensor

% Param.E = 2.563*10^5; % Young modulu
% Param.G=8.543*10^4; % shear modulu
% Param.r = 0.01; % rod radius
% Param.rho = 1.41*10^3; % Mass density
% 
% Param.g =0*[0 ; 0 ; 0 ; 9.81 ; 0; 0]; % gravity
% Param.L = 0.2; % Rod length
% Param.A = pi*Param.r^2; % Cross section area
% Param.J1 = pi*Param.r^4/2; % Polar inertia moment
% Param.J2 = pi*Param.r^4/4; % Inertia moment
% Param.J3 = pi*Param.r^4/4; % Inertia moment
% %Param.mu=1000;
% Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
% Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia
% Param.dX=Param.L/100; %spatial step
% Param.n_X=length(Param.dX:Param.dX:Param.L);
% 
% %Platform parameters
% Param.rhop=0*1.41*10^3;
% Param.a=0.02;
% Param.b=0.08;
% Param.c=0.02;
% Param.Vp=Param.a*Param.b*Param.c;
% Param.mp=Param.rhop*Param.Vp;
% %Param.Jp3=1/12*Param.mp*(Param.a^2+Param.b^2);
% %Param.Mp = diag([Param.Jp3,Param.mp,Param.mp]);%cross-sectional inertia
% Param.Jp1=1/12*Param.mp*(Param.b^2+Param.c^2);
% Param.Jp2=1/12*Param.mp*(Param.a^2+Param.c^2);
% Param.Jp3=1/12*Param.mp*(Param.a^2+Param.b^2);
% Param.Mp =diag([Param.Jp1,Param.Jp2,Param.Jp3,Param.mp,Param.mp,Param.mp]);
% %Param.Fep=[0;-Param.mp*9.81;0];
% %Param.mu=10^-3;
% 


Param.E = 2.563*10^5; % Young modulu
Param.G=8.543*10^4; % shear modulu
Param.r = 0.05; % rod radius
Param.rho = 1.41*10^3; % Mass density

Param.g =0*[0 ; 0 ; 0 ; 9.81 ; 0; 0]; % gravity
Param.L = 2; % Rod length
Param.A = pi*Param.r^2; % Cross section area
Param.J1 = pi*Param.r^4/2; % Polar inertia moment
Param.J2 = pi*Param.r^4/4; % Inertia moment
Param.J3 = pi*Param.r^4/4; % Inertia moment
%Param.mu=1000;
Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia
Param.dX=Param.L/100; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);

%Platform parametersCCC
Param.rhop=1.41*10^3;
Param.a=0.2;
Param.b=0.8;
Param.c=0.2;
Param.Vp=Param.a*Param.b*Param.c;
Param.mp=Param.rhop*Param.Vp;
Param.Jp1=1/12*Param.mp*(Param.b^2+Param.c^2);
Param.Jp2=1/12*Param.mp*(Param.a^2+Param.c^2);
Param.Jp3=1/12*Param.mp*(Param.a^2+Param.b^2);
Param.Mp =diag([Param.Jp1,Param.Jp2,Param.Jp3,Param.mp,Param.mp,Param.mp]);%cross-sectional inertia
%Param.Fep=[0;Param.mp*9.81;0];
%Param.Fep=0*[0;40;0];
%Param.mu=10^-3;
%Tendons parameters
Param.Rb= 0.25;
Param.c1=0*pi/20;
Param.c2=0*-pi/20;
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)]; % coordinate of the tendon in a cross-section
Param.n=6;  %strain modes number
Param.na=2; %number of actuated strains
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
Param.y2=0.3;
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
% 
% Param.B     = [0;0;1;0;0;0];
% Param.B_bar = [1 0 0 0 0;0 1 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;1;0;0];
% 
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
init(:,1)=zeros(2*Param.n*Param.na+Param.n_lambda+3,1);%;0;0;sqrt(Param.T_max);sqrt(Param.T_max)]];
% init(:,1)=[   0.000000000000000
%   -0.000000000000000
%   -0.000000000000019
%   -0.049999999999998
%    0.000000000000000
%   -0.000000000000000
%    0.000000000000000
%    0.000000000000000
%    0.000000000000019
%   -0.049999999999998
%    0.000000000000000
%   -0.000000000000000
%    0.000000000000000
%    4.025950985575140
%    4.025950985575142
%    0.050324387319689
%    0.000000000000002
%    0.000000000000000
%   -0.050324387319689
%   -0.000000000000002
%   -0.000000000000000
%    2.006477257677031
%    2.006477257677032
%    5.096474174802111
%    5.096474174802110];
%P_des=[0.2+Param.a/2;Param.y2/2];
plot(2+Param.a/2,Param.y2/2,'r.'); title('Workspace');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;
Ws(:,1)=[2+Param.a/2;Param.y2/2];
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
        if Neighbors(:,j)==Ws(:,l)
            t=0;
        end
    end
 if and(t==1,Neighbors(1,j)>0)
 tic
 P_des=Neighbors(:,j)
 options=optimset('fsolve');
 options.MaxFunEvals =100000000000;
 options.MaxIter =100;
 options.TolFun=1e-5;
% options.StepTolerance=10^-4;
 %options.exitflag=4;
% options.OptimalityTolerance=10^(-5);
 options.Display = 'iter';
 options.Jacobian='on';
%options.Algorithm='levenberg-marquardt';
%options.Algorithm='trust-region-dogleg';
%sing=1;
%[x,f,a,b,J1]=fsolve(@(var)Static_Gradient_Limit_actuation_V3(var,P_des,Param),init(:,i),options);
%[x,f,a,b,J1]=fsolve(@(var)Static_Gradient_Limit_actuation(var,P_des,Param),init(:,i),options);
[x,f,a,b,J1]=fsolve(@(var)Static_Gradient(var,P_des,Param),init(:,i),options);
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
q1=x(1:Param.n*Param.na);
q2=x(Param.n*Param.na+1:2*Param.n*Param.na);
qp1=x(2*Param.n*Param.na+1);
T=x(2*Param.n*Param.na+2:2*Param.n*Param.na+3)
% S1=x(2*Param.n*Param.na+1+Param.n_lambda+1:2*Param.n*Param.na+1+Param.n_lambda+2);
% S2=x(2*Param.n*Param.na+1+Param.n_lambda+3:2*Param.n*Param.na+1+Param.n_lambda+4);
    if max(abs(f))<1e-4 & T(1)>=0 & T(2)>=0 & T(1)<=100 & T(2)<=100
        hold on
        %             rx=(r1(1,end)+r2(1,end))/2;
        %             ry=(r1(2,end)+r2(2,end))/2;
%         if sing==1
        plot(P_des(1),P_des(2),'r.'); title('Workspace');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
        xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;
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
%         scatter(P_des(1),P_des(2),[],sig1,'filled')
%         colorbar
%         colormap jet
         Ws(:,k+1)=P_des;
         init(:,k+1)=x;
         k=k+1;
        %Neighbors=[Ws(:,i)+[0.005;0],Ws(:,i)+[-0.005;0],Ws(:,i)+[0;0.005],Ws(:,i)+[0;-0.005],Ws(:,i)+[0.005;0.005],Ws(:,i)+[-0.005;0.005],Ws(:,i)+[-0.005;0.005],Ws(:,i)+[-0.005;-0.005]];
        
        %Choisir le plus proche au point trouvé
        %Param.P_des=Neighbors()
    end
    end
    j=j+1;
    end
%     Param.P_des=Ws(:,i+1);
i=i+1;
    end
toc