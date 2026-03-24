%Main
%Parameters
clc
clear
tic
format long
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

Param.E = 2.563*10^5; % Young modulu
Param.G=8.543*10^4; % shear modulu
Param.r = 20; % rod radius
Param.rho = 1.41*10^3; % Mass density

Param.g =0*[0 ; 0 ; 0 ; 9.81 ; 0; 0]; % gravity
Param.L = 200; % Rod length
Param.A = pi*Param.r^2; % Cross section area
Param.J1 = pi*Param.r^4/2; % Polar inertia moment
Param.J2 = pi*Param.r^4/4; % Inertia moment
Param.J3 = pi*Param.r^4/4; % Inertia moment
%Param.mu=1000;
Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor


Param.rho = 1.41*10^3; % Mass density
Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia
Param.dX=Param.L/200; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);

%Platform parametersCCC
Param.rhop=1.41*10^3;
Param.a=20;
Param.b=130;
Param.c=20;
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
Param.Rb= 25;
Param.c1=0*pi/20;
Param.c2=0*-pi/20;
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)]; % coordinate of the tendon in a cross-section
Param.n=5;  %strain modes number
Param.na=2; %number of actuated strains
%Param.dX=Param.L/1000;%length of segments
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
  
Param.Forces_Tendons=[0
10000000
10000000
0]; 
   
Param.Ftip=0*[0;-1/sqrt(2);1/sqrt(2);0;0;0];
%Param.Forces_Tendons(1)=1;

%Initial conditions
Param.theta01=0*-pi/20;
Param.theta02=0*pi/20;
R01=[cos(Param.theta01) -sin(Param.theta01) 0; sin(Param.theta01) cos(Param.theta01) 0; 0 0 1];
R02=[cos(Param.theta02) -sin(Param.theta02) 0; sin(Param.theta02) cos(Param.theta02) 0; 0 0 1];
Param.g1=eye(4);
Param.g1(1:3,1:3)=R01;
Param.g2=eye(4);
Param.y2=80;
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
 init=zeros(2*Param.na*Param.n+3+Param.n_lambda,1);
%   init=[  -4.160536578227304
%    0.101607374549173
%    9.510238601500443
%    0.064295054247374
%   -1.292727856487818
%   -0.010023823641902
%   -0.238207191358419
%    0.000207440110509
%    0.003958114599368
%    0.000044764567038
%   -0.003398132182113
%   -0.000066147726476
%   -4.163965732925931
%   -0.278837512751037
%   -5.225948137397218
%   -0.038878222678497
%   -0.729300699526419
%    0.002096593901251
%    0.054498482766487
%   -0.000122486708428
%   -0.002543226202195
%   -0.000067690662247
%   -0.001063089680145
%   -0.000008723733428
%   -0.832315043908110
%    0.167109380905745
%   -0.037512122069738
%    0.819042013869811
%  -16.391293630860382
%    7.071533512671114
%   -0.168506176778937
%   16.391293630860382
%   -7.071533512671114];
 init(2*Param.n*Param.na+2:2*Param.n*Param.na+3)=[200+Param.a/2;Param.y2/2];
 options=optimset('fsolve');
 %options.MaxFunEvals =1;
 options.MaxIter = 10000;
 options.TolFun=1e-5;
 %options.StepTolerance=10^-5;
 %options.exitflag=4;
%options.OptimalityTolerance=10^(-5);
 options.Display = 'iter';
options.Jacobian='on';
%options.Algorithm='levenberg-marquardt';
%options.Algorithm='trust-region-dogleg';
[x,f,a,b,J]=fsolve(@(var)Static_gradient(var,Param),init,options);
%options.Jacobian='off';
%[x,f,a,b,J1]=fsolve(@(var)Static(var,Param),init,options);
%[x,f,a,b,J2]=fsolve(@(var)Static_V2(var,Param),init,options);
%J
%J1
%J-J1
q1=x(1:Param.n*Param.na);
q2=x(Param.n*Param.na+1:2*Param.n*Param.na);
qp=x(2*Param.n*Param.na+1:2*Param.n*Param.na+3);
norm(f)
toc

r1=zeros(3,Param.n_X+1);
r2=zeros(3,Param.n_X+1);
r2(2,1)=Param.y2;
Q1=zeros(4,Param.n_X+1);
Q2=zeros(4,Param.n_X+1);
Q1(:,1)=[cos(Param.theta01/2);0;0;sin(Param.theta01/2)];
Q2(:,1)=[cos(Param.theta02/2);0;0;sin(Param.theta02/2)];
for j=1:Param.n_X
  %   g=Geometric_model((j-1)*Param.dX,q,Param);
  %   r(:,j)=g(1:3,4);
        phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
        xia1=Param.xi_a0 + phi*q1;
        xi1=Param.B*xia1+Param.B_bar*Param.xi_c;
        xia2=Param.xi_a0 + phi*q2;
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
        r1(:,j+1)=r1(:,j) + r_X1*Param.dX;
        R1=eye(3) + 2/(Q1(:,j+1)'*Q1(:,j+1)) * [-Q1(3,j+1)^2-Q1(4,j+1)^2, Q1(2,j+1)*Q1(3,j+1)-Q1(4,j+1)*Q1(1,j+1),Q1(2,j+1)*Q1(4,j+1) + Q1(3,j+1)*Q1(1,j+1) ; 
                                 Q1(2,j+1)*Q1(3,j+1)+Q1(4,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(4,j+1)^2,Q1(3,j+1)*Q1(4,j+1) - Q1(2,j+1)*Q1(1,j+1) ; 
                                 Q1(2,j+1)*Q1(4,j+1)-Q1(3,j+1)*Q1(1,j+1), Q1(3,j+1)*Q1(4,j+1) + Q1(2,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(3,j+1)^2];

         
        R2=eye(3) + 2/(Q2(:,j)'*Q2(:,j)) * [-Q2(3,j)^2-Q2(4,j)^2, Q2(2,j)*Q2(3,j)-Q2(4,j)*Q2(1,j),Q2(2,j)*Q2(4,j) + Q2(3,j)*Q2(1,j) ; 
                                 Q2(2,j)*Q2(3,j)+Q2(4,j)*Q2(1,j), -Q2(2,j)^2-Q2(4,j)^2,Q2(3,j)*Q2(4,j) - Q2(2,j)*Q2(1,j) ; 
                                 Q2(2,j)*Q2(4,j)-Q2(3,j)*Q2(1,j), Q2(3,j)*Q2(4,j) + Q2(2,j)*Q2(1,j), -Q2(2,j)^2-Q2(3,j)^2];
        Q_X2 = [ 0, -K2(1), -K2(2), -K2(3);
                 K2(1), 0, K2(3), -K2(2);
                K2(2), -K2(3), 0, K2(1);
                K2(3), K2(2), -K2(1), 0 ] * Q2(:,j)/2;
        r_X2 = R2*Gamma2;
        Q2(:,j+1)=Q2(:,j) + Q_X2*Param.dX;
        r2(:,j+1)=r2(:,j) + r_X2*Param.dX;
        R2=eye(3) + 2/(Q2(:,j+1)'*Q2(:,j+1)) * [-Q2(3,j+1)^2-Q2(4,j+1)^2, Q2(2,j+1)*Q2(3,j+1)-Q2(4,j+1)*Q2(1,j+1),Q2(2,j+1)*Q2(4,j+1) + Q2(3,j+1)*Q2(1,j+1) ; 
                                 Q2(2,j+1)*Q2(3,j+1)+Q2(4,j+1)*Q2(1,j+1), -Q2(2,j+1)^2-Q2(4,j+1)^2,Q2(3,j+1)*Q2(4,j+1) - Q2(2,j+1)*Q2(1,j+1) ; 
                                 Q2(2,j+1)*Q2(4,j+1)-Q2(3,j+1)*Q2(1,j+1), Q2(3,j+1)*Q2(4,j+1) + Q2(2,j+1)*Q2(1,j+1), -Q2(2,j+1)^2-Q2(3,j+1)^2];
end  
g1=[R1 r1(:,end);0 0 0 1];
g2=[R2 r2(:,end);0 0 0 1];
Rp=[cos(qp(1)) -sin(qp(1)); sin(qp(1)) cos(qp(1))];
         hold on
    l_rec=0.02;
    plot(r1(1,:),r1(2,:),'b','LineWidth', 4); title('Planar TACR');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    plot(r2(1,:),r2(2,:),'b','LineWidth', 4);
        plot(qp(2),qp(3),'g*','LineWidth', 4);
%       xrec=[qp(2); qp(2)+2*Param.a/2*cos(qp(1))];
%       yrec=[qp(3); qp(3)+2*Param.a/2*sin(qp(1))];
%       plot(xrec,yrec,'r','LineWidth', 4);
     xrec=[qp(2)+(Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1))); qp(2)+(Param.a/2*cos(qp(1))+Param.b/2*sin(qp(1))) ;qp(2)-(Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1)));qp(2)+(-Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1)));qp(2)+(Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1)))];
     yrec=[qp(3)+(+Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1))) ;qp(3)+(Param.a/2*sin(qp(1))-Param.b/2*cos(qp(1))) ;qp(3)-(Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1)));qp(3)+(-Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1)));qp(3)+(Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1)))];
     plot(xrec,yrec,'r','LineWidth', 4);
   % quiver3(r(1,end),r(2,end),r(3,end),r(1,end)+Ft(1),r(2,end)+Ft(2),r(3,end)+Ft(3));