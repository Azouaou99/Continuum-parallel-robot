%Main
%Parameters
clc
clear
format long
%Param.E = 3.5*10^6; % Young modulu
Param.E = 3*10^9; % Young modulu
Param.G=8.543*10^4; % shear modulu
%Param.r = 0.02; % rod radius
Param.rho = 1240*0.2; % Mass density
mp=0;
dis=0.04;
u=1.5e-4;
%Platform parameters
Param.a=0.02;
Param.b=0.1;
Param.c=0.03;
Param.rhop=1.41*10^3;
Param.Vp=Param.a*Param.b*Param.c;
Param.mp=mp;
Param.Jp1=1/12*Param.mp*(Param.b^2+Param.c^2);
Param.Jp2=1/12*Param.mp*(Param.a^2+Param.c^2);
Param.Jp3=1/12*Param.mp*(Param.a^2+Param.b^2);
Param.Mp =diag([Param.Jp3,Param.mp,Param.mp]);%cross-sectional inertia
Param.Mp2 =diag([Param.Jp1,Param.Jp2,Param.Jp3,Param.mp,Param.mp,Param.mp]);%cross-sectional inertia
Param.F_p=[0;Param.mp*9.81;0];

Param.g =[0 ; 0 ; 0 ; 9.81 ; 0; 0]; % gravity
Param.L = 0.2; % Rod length
Param.d=0.001;
Param.A = Param.c*Param.d; % Cross section area
Param.J1 = Param.c*Param.d*(Param.c^2+Param.d^2)/12; % Polar inertia moment
Param.J2 = Param.c^3*Param.d/12; % Inertia moment
Param.J3 = Param.c*Param.d^3/12; % Inertia moment
%Param.mu=1000;
%u=10000*Param.J3/Param.A
Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.J3/u, Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia
Param.dX=Param.L/100; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);
%Param.mu=0.4;
%Param.D = Param.mu*diag([Param.J1 , 3*Param.J2 , 3*Param.J3 , 3*Param.A , Param.A , Param.A]);%Damping matrix
Param.D=10^(-4)*eye(6);%Damping matrix
Param.duree=30;%simulation time
Param.dt=0.01; %time step
Param.dX=Param.L/100; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);
Param.n_t=length(Param.dt:Param.dt:Param.duree);
Param.n_lambda=6;



%Tendons parameters
Param.Rb= 0.015; % Distance between a tendon and the backbone
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)];% coordinate of the tendon in a cross-section
%Param.Tendons_list = [Param.Tendon_coordinate(0),Param.Tendon_coordinate(pi/2),Param.Tendon_coordinate(pi),Param.Tendon_coordinate(3*pi/2)]; % 4 tendons separated with an angle of pi/2
Param.n=3;  %strain modes number
Param.na=2; %number of actuated strains

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
%Param.Forces_Tendons = zeros(4,1); % tension in each tendon
   
   

%   
Param.Forces_Tendons=[0
   %0.035
   12.1
  0
  0
];   
% 
% Param.Forces_Tendons=[0
% 0
% 0
% 0]; 
% %Param.Ftip=[0;0;0;0;0;2];

%Initial conditions
Param.c1=0*pi/400;
Param.c2=0*-pi/400;
Param.theta01=0*pi/8;
Param.theta02=0*-pi/8;
R01=[cos(Param.theta01) -sin(Param.theta01) 0; sin(Param.theta01) cos(Param.theta01) 0; 0 0 1];
R02=[cos(Param.theta02) -sin(Param.theta02) 0; sin(Param.theta02) cos(Param.theta02) 0; 0 0 1];
Param.g1=eye(4);
Param.g1(1:3,1:3)=R01;
Param.g2=eye(4);
Param.y2=dis;
Param.g2(1:3,1:3)=R02;
Param.g2(2,4)=Param.y2;
Param.CM1=[Param.a/2;Param.y2/2];
Param.CM2=[Param.a/2;-Param.y2/2];

% Param.B =[0 0;1 0;0 0;0 0;0 0;0 1];
% Param.B_bar = [1 0 0 0;0 0 0 0; 0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;0;1];

Param.B     = [0 0;0 0;1 0;0 1;0 0;0 0];%Selection matrix 
Param.B_bar =[1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0;0 0 1 0;0 0 0 1];
Param.xi_0=[0;0;0;1;0;0];%reference configuration
Param.xi_a0=Param.B'*Param.xi_0; Param.xi_c=[0;0;0;0];

% Param.B     = [0;0;1;0;0;0];%Selection matrix
% Param.B_bar = [1 0 0 0 0;0 1 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
% Param.xi_0=[0;0;0;1;0;0];%reference configuration
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;1;0;0];

% Param.B     = [0 0 0 ;0 0 0;1 0 0;0 1 0;0 0 1;0 0 0];
% Param.B_bar = [1 0 0;0 1 0;0 0 0;0 0 0;0 0 0;0 0 1];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;0];

Param.Ha=Param.B'*Param.H*Param.B;%the matrix of the reduced Hooke coefficients
Param.Da=Param.B'*Param.D*Param.B;%reduced damping matrix

  
  
  

  


 
%Initial conditions
%  Param.q10=[ -1.504989513688684
%   -0.000000000000000
%   18.155268566737533
%   -0.000132577106367
%   -0.000000000000000
%    0.000001406257921];
% % 
%  Param.q20=[-1.504989513688684
%   -0.000000000000001
%   -2.332743627614550
%    0.000069442963581
%    0.000000000000000
%   -0.000000409061144];
% % 
% Param.qp0=[  -0.300997902737737
%    0.201877600262403
%   -0.015130306863276];
  
 
 
 

x0=zeros(2*Param.n*Param.na+3+Param.n_lambda,1);
% x0=[  -2.601050783392311
%   -0.020688135768422
%   10.991590589140579
%   -0.054828594182420
%   -0.000023679751450
%    0.000553266793621
%   -2.601050783372423
%   -0.018707423954267
%   -3.403709568592795
%    0.015417855323013
%   -0.000025492298079
%   -0.000282775820591
%   -0.520210156675189
%    0.193753732675140
%   -0.034094785256798
%    0.114677292096575
%   -4.492736552504298
%    1.192328317622690
%   -0.024855783696523
%    4.492736552504298
%   -1.192328317622690];

Param.q10=x0(1:Param.n*Param.na);
Param.q20=x0(Param.n*Param.na+1:2*Param.n*Param.na);
%Param.qp0=x0(2*Param.n*Param.na+1:2*Param.n*Param.na+3);
Param.qp0=[0;0.2+Param.a/2;Param.y2/2];
lambda=zeros(Param.n_lambda,Param.n_t+1);
lambda(:,1)=x0(2*Param.n*Param.na+4 : 2*Param.n*Param.na+3+Param.n_lambda);
Param.q0=[Param.q10;Param.q20;Param.qp0];
Param.qdp=zeros(3,1);
Param.qd0=zeros(Param.n*Param.na,1);
%initialization of the coefficients
q1=zeros(Param.na*Param.n,Param.n_t+1);
q2=zeros(Param.na*Param.n,Param.n_t+1);
qp=zeros(3,Param.n_t+1);
q=zeros(2*Param.na*Param.n+3,Param.n_t+1);

q1(:,1)=Param.q10;
q2(:,1)=Param.q20;
qp(:,1)=Param.qp0;

q(:,1)=Param.q0;
qd1=zeros(Param.na*Param.n,Param.n_t+1);
qd2=zeros(Param.na*Param.n,Param.n_t+1);
qdp=zeros(3,Param.n_t+1);
qd1(:,1)=Param.qd0;
qd2(:,1)=Param.qd0;
qdp(:,1)=Param.qdp;
qd=zeros(2*Param.na*Param.n+3,Param.n_t+1);


qdd1=zeros(Param.na*Param.n,Param.n_t+1);
qdd2=zeros(Param.na*Param.n,Param.n_t+1);
qddp=zeros(3,Param.n_t+1);


%r=zeros(3,Param.n_X+1,Param.n_t+1);
%Q=zeros(4,Param.n_X+1);
%Q(:,1)=[1;0;0;0];

r1=zeros(3,Param.n_X+1,Param.n_t+1);
r2=zeros(3,Param.n_X+1,Param.n_t+1);
r2(2,1,:)=Param.y2;
% for j=1:Param.n_X
%     phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
%     xia=Param.xi_a0 + phi*q(:,1);
%     xi=Param.B*xia+Param.B_bar*Param.xi_c;
%     K=xi(1:3); %angular strain
%     Gamma=xi(4:6);%Linear strain
%     R=eye(3) + 2/(Q(:,j)'*Q(:,j)) * [-Q(3,j)^2-Q(4,j)^2, Q(2,j)*Q(3,j)-Q(4,j)*Q(1,j),Q(2,j)*Q(4,j) + Q(3,j)*Q(1,j) ;
%         Q(2,j)*Q(3,j)+Q(4,j)*Q(1,j), -Q(2,j)^2-Q(4,j)^2,Q(3,j)*Q(4,j) - Q(2,j)*Q(1,j) ;
%         Q(2,j)*Q(4,j)-Q(3,j)*Q(1,j), Q(3,j)*Q(4,j) + Q(2,j)*Q(1,j), -Q(2,j)^2-Q(3,j)^2];
%     Q_X = [ 0, -K(1), -K(2), -K(3);
%         K(1), 0, K(3), -K(2);
%         K(2), -K(3), 0, K(1);
%         K(3), K(2), -K(1), 0 ] * Q(:,j)/2;
%     r_X = R*Gamma;
%     Q(:,j+1)=Q(:,j) + Q_X*Param.dX;
%     r(:,j+1,1)=r(:,j,1) + r_X*Param.dX;
% end
Q1=zeros(4,Param.n_X+1);
Q2=zeros(4,Param.n_X+1);
Q1(:,1)=[cos(Param.theta01/2);0;0;sin(Param.theta01/2)];
Q2(:,1)=[cos(Param.theta02/2);0;0;sin(Param.theta02/2)];
for j=1:Param.n_X
    phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
    xia1=Param.xi_a0 + phi*q1(:,1);
    xi1=Param.B*xia1+Param.B_bar*Param.xi_c;
    xia2=Param.xi_a0 + phi*q2(:,1);
    xi2=Param.B*xia2+Param.B_bar*Param.xi_c;
    %xi(4:6)=[1;0;0];
    K1=xi1(1:3); %angular strain
    Gamma1=xi1(4:6);%Linear strain
    K2=xi2(1:3); %angular strain
    Gamma2=xi2(4:6);%Linear strain

    R1=eye(3) + 2/(Q1(:,j)'*Q1(:,j)) * [-Q1(3,j)^2-Q1(4,j)^2, Q1(2,j)*Q1(3,j)-Q1(4,j)*Q1(1,j),Q1(2,j)*Q1(4,j) + Q1(3,j)*Q1(1,j) ;
        Q1(2,j)*Q1(3,j)+Q1(4,j)*Q1(1,j), -Q1(2,j)^2-Q1(4,j)^2,Q1(3,j)*Q1(4,j) - Q1(2,j)*Q1(1,j) ;
        Q1(2,j)*Q1(4,j)-Q1(3,j)*Q1(1,j), Q1(3,j)*Q1(4,j) + Q1(2,j)*Q1(1,j), -Q1(2,j)^2-Q1(3,j)^2];
    Q_X1 = [ 0, -K1(1), -K1(2), -K1(3);
        K1(1), 0, K1(3), -K1(2);
        K1(2), -K1(3), 0, K1(1);
        K1(3), K1(2), -K1(1), 0 ] * Q1(:,j)/2;
    r_X1 = R1*Gamma1;
    Q1(:,j+1)=Q1(:,j) + Q_X1*Param.dX;
    r1(:,j+1,1)=r1(:,j,1) + r_X1*Param.dX;
    %  R1=eye(3) + 2/(Q1(:,j+1)'*Q1(:,j+1)) * [-Q1(3,j+1)^2-Q1(4,j+1)^2, Q1(2,j+1)*Q1(3,j+1)-Q1(4,j+1)*Q1(1,j+1),Q1(2,j+1)*Q1(4,j+1) + Q1(3,j+1)*Q1(1,j+1) ;
    %                           Q1(2,j+1)*Q1(3,j+1)+Q1(4,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(4,j+1)^2,Q1(3,j+1)*Q1(4,j+1) - Q1(2,j+1)*Q1(1,j+1) ;
    %                           Q1(2,j+1)*Q1(4,j+1)-Q1(3,j+1)*Q1(1,j+1), Q1(3,j+1)*Q1(4,j+1) + Q1(2,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(3,j+1)^2];


    R2=eye(3) + 2/(Q2(:,j)'*Q2(:,j)) * [-Q2(3,j)^2-Q2(4,j)^2, Q2(2,j)*Q2(3,j)-Q2(4,j)*Q2(1,j),Q2(2,j)*Q2(4,j) + Q2(3,j)*Q2(1,j) ;
        Q2(2,j)*Q2(3,j)+Q2(4,j)*Q2(1,j), -Q2(2,j)^2-Q2(4,j)^2,Q2(3,j)*Q2(4,j) - Q2(2,j)*Q2(1,j) ;
        Q2(2,j)*Q2(4,j)-Q2(3,j)*Q2(1,j), Q2(3,j)*Q2(4,j) + Q2(2,j)*Q2(1,j), -Q2(2,j)^2-Q2(3,j)^2];
    Q_X2 = [ 0, -K2(1), -K2(2), -K2(3);
        K2(1), 0, K2(3), -K2(2);
        K2(2), -K2(3), 0, K2(1);
        K2(3), K2(2), -K2(1), 0 ] * Q2(:,j)/2;
    r_X2 = R2*Gamma2;
    Q2(:,j+1)=Q2(:,j) + Q_X2*Param.dX;
    r2(:,j+1,1)=r2(:,j,1) + r_X2*Param.dX;
    %R2=eye(3) + 2/(Q2(:,j+1)'*Q2(:,j+1)) * [-Q2(3,j+1)^2-Q2(4,j+1)^2, Q2(2,j+1)*Q2(3,j+1)-Q2(4,j+1)*Q2(1,j+1),Q2(2,j+1)*Q2(4,j+1) + Q2(3,j+1)*Q2(1,j+1) ;
    %                         Q2(2,j+1)*Q2(3,j+1)+Q2(4,j+1)*Q2(1,j+1), -Q2(2,j+1)^2-Q2(4,j+1)^2,Q2(3,j+1)*Q2(4,j+1) - Q2(2,j+1)*Q2(1,j+1) ;
    %                         Q2(2,j+1)*Q2(4,j+1)-Q2(3,j+1)*Q2(1,j+1), Q2(3,j+1)*Q2(4,j+1) + Q2(2,j+1)*Q2(1,j+1), -Q2(2,j+1)^2-Q2(3,j+1)^2];
end

Y=0:Param.dX:Param.L;
SK=zeros(Param.n*Param.na,Param.n*Param.na);
SD=zeros(Param.n*Param.na,Param.n*Param.na);
for k=1:Param.n_X+1
    phival=Phi(Param.na,Param.n,Y(k),Param.L);
    fK=phival'*Param.Ha*phival;
    fD=phival'*Param.Da*phival;
    if k==1
        fK0=fK;
        fD0=fD;
    end
    if k==Param.n_X+1
        fKn=fK;
        fDn=fD;
    end
    if and(k~=Param.n_X+1,k~=1)
        SK=SK+fK;
        SD=SD+fD;
    end
end
Param.Keps= Param.L/(Param.n_X)*((fK0+fKn)/2 + SK);%Stiffness matrix on the modal space
Param.Deps= Param.L/(Param.n_X)*((fD0+fDn)/2 + SD);%Damping matrix on the modal space

Param.Y=0:Param.dX:Param.L;
l=1;
% T=15;
   for i=1:Param.n_t %Time updating
%        %T=T+0.001;
%        %T=T+0.05;
%        if i>50
%    Param.Forces_Tendons=[0
%  T
%  0
%  0];
%       end
%       if i>Param.n_t/2 &i<Param.n_t/2+10
%  Param.Ftendons=[1;1;0;1];
%       else 
%        Param.Ftendons=[1;0;0;1];   
%      end
    disp (i*Param.dt)
    init=[q1(:,i);q2(:,i);qp(:,i);lambda(:,i)];
   % options = optimoptions("fsolve");


%      options=optimset('fsolve');
%      options.MaxFunEvals =1000000000000;
%      options.MaxIter = 50000000000;
%      options.TolFun=1e-7;
%      options.Jacobian='on';
% %     % options.StepTolerance=10^-4;
% %     %options.exitflag=4;
% %     %options.OptimalityTolerance=10^(-2);
%      options.Display = 'iter';
%     %options.Algorithm='levenberg-marquardt';
        % 
        % options  =  optimoptions( 'fsolve', ...
        %                   'FiniteDifferenceStepSize', 1e-7, ...
        %                   'OptimalityTolerance', 1e-7, ...
        %                   'StepTolerance', 1e-7, ...
        %                   'FunctionTolerance', 1e-7, ...
        %                   'SpecifyObjectiveGradient', true, ...
        %                   'Display', 'off');
    %x=fsolve(@(var)Dynamic_grad(var,q(:,i),qd(:,i),i,Param),init,options);%2nd derivative of q
    tic    
x=Newton_r(init,q(:,i),qd(:,i),1e-4,Param);
toc
    %Param.Forces_Tendons(1,1)=Param.Forces_Tendons(1,1)+0.01;
   % Param.Forces_Tendons(1,2)=Param.Forces_Tendons(1,2)+0.01;
    %Param.Forces_Tendons(2,3)=Param.Forces_Tendons(2,3)+0.01;
   % Param.Forces_Tendons(2,4)=Param.Forces_Tendons(2,4)+0.01;
%     qd1(:,i+1)=x(1:Param.n*Param.na);
%     qd2(:,i+1)=x(Param.n*Param.na+1:2*Param.n*Param.na);
%     qd(:,i+1)=x(1:2*Param.n*Param.na);
    q1(:,i+1)=x(1:Param.n*Param.na);
    q2(:,i+1)=x(Param.n*Param.na+1:2*Param.n*Param.na);
    qp(:,i+1)=x(2*Param.n*Param.na+1:2*Param.n*Param.na+3);
    q(:,i+1)=x(1:2*Param.n*Param.na+3);
    lambda(:,i+1)=x(2*Param.n*Param.na+4 : 2*Param.n*Param.na+3+Param.n_lambda);

    qd1(:,i+1)=(q1(:,i+1)-q1(:,i))/Param.dt;%x(2*Param.n*Param.na+1:3*Param.n*Param.na);%time intgeration
    qd2(:,i+1)=(q2(:,i+1)-q2(:,i))/Param.dt;%x(3*Param.n*Param.na+1:4*Param.n*Param.na);
    qdp(:,i+1)=(qp(:,i+1)-qp(:,i))/Param.dt;%x(3*Param.n*Param.na+1:4*Param.n*Param.na);
    qd(:,i+1)=(q(:,i+1)-q(:,i))/Param.dt;

    Q1(:,1)=[cos(Param.theta01/2);0;0;sin(Param.theta01/2)];
    Q2(:,1)=[cos(Param.theta02/2);0;0;sin(Param.theta02/2)];
    for j=1:Param.n_X
        phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
        xia1=Param.xi_a0 + phi*q1(:,i+1);
        xi1=Param.B*xia1+Param.B_bar*Param.xi_c;
        xia2=Param.xi_a0 + phi*q2(:,i+1);
        xi2=Param.B*xia2+Param.B_bar*Param.xi_c;
        %xi(4:6)=[1;0;0];
        K1=xi1(1:3); %angular strain
        Gamma1=xi1(4:6);%Linear strain
        K2=xi2(1:3); %angular strain
        Gamma2=xi2(4:6);%Linear strain

        R1=eye(3) + 2/(Q1(:,j)'*Q1(:,j)) * [-Q1(3,j)^2-Q1(4,j)^2, Q1(2,j)*Q1(3,j)-Q1(4,j)*Q1(1,j),Q1(2,j)*Q1(4,j) + Q1(3,j)*Q1(1,j) ;
            Q1(2,j)*Q1(3,j)+Q1(4,j)*Q1(1,j), -Q1(2,j)^2-Q1(4,j)^2,Q1(3,j)*Q1(4,j) - Q1(2,j)*Q1(1,j) ;
            Q1(2,j)*Q1(4,j)-Q1(3,j)*Q1(1,j), Q1(3,j)*Q1(4,j) + Q1(2,j)*Q1(1,j), -Q1(2,j)^2-Q1(3,j)^2];
        Q_X1 = [ 0, -K1(1), -K1(2), -K1(3);
            K1(1), 0, K1(3), -K1(2);
            K1(2), -K1(3), 0, K1(1);
            K1(3), K1(2), -K1(1), 0 ] * Q1(:,j)/2;
        r_X1 = R1*Gamma1;
        Q1(:,j+1)=Q1(:,j) + Q_X1*Param.dX;
        r1(:,j+1,i+1)=r1(:,j,i+1) + r_X1*Param.dX;
        %  R1=eye(3) + 2/(Q1(:,j+1)'*Q1(:,j+1)) * [-Q1(3,j+1)^2-Q1(4,j+1)^2, Q1(2,j+1)*Q1(3,j+1)-Q1(4,j+1)*Q1(1,j+1),Q1(2,j+1)*Q1(4,j+1) + Q1(3,j+1)*Q1(1,j+1) ;
        %                           Q1(2,j+1)*Q1(3,j+1)+Q1(4,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(4,j+1)^2,Q1(3,j+1)*Q1(4,j+1) - Q1(2,j+1)*Q1(1,j+1) ;
        %                           Q1(2,j+1)*Q1(4,j+1)-Q1(3,j+1)*Q1(1,j+1), Q1(3,j+1)*Q1(4,j+1) + Q1(2,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(3,j+1)^2];


        R2=eye(3) + 2/(Q2(:,j)'*Q2(:,j)) * [-Q2(3,j)^2-Q2(4,j)^2, Q2(2,j)*Q2(3,j)-Q2(4,j)*Q2(1,j),Q2(2,j)*Q2(4,j) + Q2(3,j)*Q2(1,j) ;
            Q2(2,j)*Q2(3,j)+Q2(4,j)*Q2(1,j), -Q2(2,j)^2-Q2(4,j)^2,Q2(3,j)*Q2(4,j) - Q2(2,j)*Q2(1,j) ;
            Q2(2,j)*Q2(4,j)-Q2(3,j)*Q2(1,j), Q2(3,j)*Q2(4,j) + Q2(2,j)*Q2(1,j), -Q2(2,j)^2-Q2(3,j)^2];
        Q_X2 = [ 0, -K2(1), -K2(2), -K2(3);
            K2(1), 0, K2(3), -K2(2);
            K2(2), -K2(3), 0, K2(1);
            K2(3), K2(2), -K2(1), 0 ] * Q2(:,j)/2;
        r_X2 = R2*Gamma2;
        Q2(:,j+1)=Q2(:,j) + Q_X2*Param.dX;
        r2(:,j+1,i+1)=r2(:,j,i+1) + r_X2*Param.dX;
        %R2=eye(3) + 2/(Q2(:,j+1)'*Q2(:,j+1)) * [-Q2(3,j+1)^2-Q2(4,j+1)^2, Q2(2,j+1)*Q2(3,j+1)-Q2(4,j+1)*Q2(1,j+1),Q2(2,j+1)*Q2(4,j+1) + Q2(3,j+1)*Q2(1,j+1) ;
        %                         Q2(2,j+1)*Q2(3,j+1)+Q2(4,j+1)*Q2(1,j+1), -Q2(2,j+1)^2-Q2(4,j+1)^2,Q2(3,j+1)*Q2(4,j+1) - Q2(2,j+1)*Q2(1,j+1) ;
        %                         Q2(2,j+1)*Q2(4,j+1)-Q2(3,j+1)*Q2(1,j+1), Q2(3,j+1)*Q2(4,j+1) + Q2(2,j+1)*Q2(1,j+1), -Q2(2,j+1)^2-Q2(3,j+1)^2];
    end


% l_rec=0.02;
% for j=1:Param.n_t+1
%    % j=50
%     plot3(r1(1,:,j),r1(2,:,j),r1(3,:,j),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
%     xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
%     hold on
%     plot3(r2(1,:,j),r2(2,:,j),r2(3,:,j),'b','LineWidth', 4); title('Cosserat rod');drawnow
%     %xrec=[r1(1,end,j);r1(1,end,j)+R1(1,1)*l_rec ; r2(1,end)+R2(1,1)*l_rec ; r2(1,end) ; r1(1,end)];
%     % yrec=[r1(2,end);r1(2,end)+R1(2,1)*l_rec;r2(2,end)+R2(2,1)*l_rec;r2(2,end);r1(2,end)];
%     % plot(xrec,yrec,'r','LineWidth', 4)
%     hold off
% end

l_rec=0.02;
%for j=1:Param.n_t+1
   % j=50
    plot(r1(1,:,i),r1(2,:,i),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
    hold on
    plot(r2(1,:,i),r2(2,:,i),'b','LineWidth', 4); title('Cosserat rod');drawnow
    plot(qp(2,i),qp(3,i),'g*','LineWidth', 4); title('Cosserat rod');drawnow
     xrec=[qp(2,i)+(Param.a/2*cos(qp(1,i))-Param.b/2*sin(qp(1,i))); qp(2,i)+(Param.a/2*cos(qp(1,i))+Param.b/2*sin(qp(1,i))) ;qp(2,i)-(Param.a/2*cos(qp(1,i))-Param.b/2*sin(qp(1,i)));qp(2,i)+(-Param.a/2*cos(qp(1,i))-Param.b/2*sin(qp(1,i)));qp(2,i)+(Param.a/2*cos(qp(1,i))-Param.b/2*sin(qp(1,i)))];
     yrec=[qp(3,i)+(+Param.a/2*sin(qp(1,i))+Param.b/2*cos(qp(1,i))) ;qp(3,i)+(Param.a/2*sin(qp(1,i))-Param.b/2*cos(qp(1,i))) ;qp(3,i)-(Param.a/2*sin(qp(1,i))+Param.b/2*cos(qp(1,i)));qp(3,i)+(-Param.a/2*sin(qp(1,i))+Param.b/2*cos(qp(1,i)));qp(3,i)+(Param.a/2*sin(qp(1,i))+Param.b/2*cos(qp(1,i)))];
     plot(xrec,yrec,'r','LineWidth', 4);drawnow
    %xrec=[r1(1,end,j);r1(1,end,j)+R1(1,1)*l_rec ; r2(1,end)+R2(1,1)*l_rec ; r2(1,end) ; r1(1,end)];
    % yrec=[r1(2,end);r1(2,end)+R1(2,1)*l_rec;r2(2,end)+R2(2,1)*l_rec;r2(2,end);r1(2,end)];
    % plot(xrec,yrec,'r','LineWidth', 4)
    hold off
%end
%  j=1;
%  plot(r1(1,:,j),r1(2,:,j),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
%     xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
%     hold on
%     plot(r2(1,:,j),r2(2,:,j),'b','LineWidth', 4); title('Cosserat rod');drawnow
%     plot(qp(2,j),qp(3,j),'g*','LineWidth', 4); title('Cosserat rod');;drawnow
%      xrec=[qp(2,j)+(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j))); qp(2,j)+(Param.a/2*cos(qp(1,j))+Param.b/2*sin(qp(1,j))) ;qp(2,j)-(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)));qp(2,j)+(-Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)));qp(2,j)+(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)))];
%      yrec=[qp(3,j)+(+Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j))) ;qp(3,j)+(Param.a/2*sin(qp(1,j))-Param.b/2*cos(qp(1,j))) ;qp(3,j)-(Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)));qp(3,j)+(-Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)));qp(3,j)+(Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)))];
%      plot(xrec,yrec,'r','LineWidth', 4);drawnow
% im={}; 
% [im{j},map]=frame2im(getframe);
% for j=2:Param.n_t+1
%     figure(1)
%     clf
%     disp(j);
%    plot(r1(1,:,j),r1(2,:,j),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
%     xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
%     hold on
%     plot(r2(1,:,j),r2(2,:,j),'b','LineWidth', 4); title('Cosserat rod');drawnow
%     plot(qp(2,j),qp(3,j),'g*','LineWidth', 4); title('Cosserat rod');;drawnow
%      xrec=[qp(2,j)+(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j))); qp(2,j)+(Param.a/2*cos(qp(1,j))+Param.b/2*sin(qp(1,j))) ;qp(2,j)-(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)));qp(2,j)+(-Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)));qp(2,j)+(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)))];
%      yrec=[qp(3,j)+(+Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j))) ;qp(3,j)+(Param.a/2*sin(qp(1,j))-Param.b/2*cos(qp(1,j))) ;qp(3,j)-(Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)));qp(3,j)+(-Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)));qp(3,j)+(Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)))];
%      plot(xrec,yrec,'r','LineWidth', 4);drawnow
%     %xrec=[r1(1,end,j);r1(1,end,j)+R1(1,1)*l_rec ; r2(1,end)+R2(1,1)*l_rec ; r2(1,end) ; r1(1,end)];
%     % yrec=[r1(2,end);r1(2,end)+R1(2,1)*l_rec;r2(2,end)+R2(2,1)*l_rec;r2(2,end);r1(2,end)];
% %     % plot(xrec,yrec,'r','LineWidth', 4)
%     [im{j},map]=frame2im(getframe);
% hold off;
%  pause(0.001)
% end
% [temp,map]=rgb2ind(im{1},Param.n_t+1);
% for j=1:Param.n_t+1
%     gifim(:,:,1,j)=rgb2ind(im{j},map);
% end
% imwrite(gifim,map,'PTACR5.gif');
end

l_ana=zeros(Param.n_t+1,1);
m_ana=zeros(Param.n_t+1,1);
n_ana=zeros(Param.n_t+1,1);
t=(0:Param.dt:Param.duree)';
for j=1:Param.n_t+1
    l_ana(j,1)=(r1(1,end,j)+r2(1,end,j))/2;
    m_ana(j,1)=(r1(2,end,j)+r2(2,end,j))/2;
    n_ana(j,1)=(r1(3,end,j)+r2(3,end,j))/2;
end
figure; hold on
grid
plot(t',l_ana,'r')
plot(t',m_ana,'b')
plot(t',n_ana,'g')
ylabel('$x$, $y$, $z$ [m] ','interpreter', 'latex');
xlabel('$t$ [s]','interpreter', 'latex');
legend({'$x_{eff}$', '$y_{eff}$', '$z_{eff}$'}, 'interpreter', 'latex')
title('End effector positions')

% figure; hold on
% grid
% plot(t',qp(2,:),'r')
% plot(t',qp(3,:),'b')
% %plot(t',n_ana,'g')
% ylabel('$x$, $y$ [m] ','interpreter', 'latex');
% xlabel('$t$ [s]','interpreter', 'latex');
% legend({'$x_{p}$', '$y_{p}$'}, 'interpreter', 'latex')
% title('End effector positions')


