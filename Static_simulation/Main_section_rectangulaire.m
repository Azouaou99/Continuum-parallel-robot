%Main
%Parameters
clc
clear
tic
format long
%  Param.E = 10^10; % Young modulu
% Param.r = 0.005; % rod radius
% Param.g =0*[0 ; 0 ; 0 ; 0 ; 0; -9.81]; % gravity
% Param.L = 1; % Rod length
% Param.A =  pi*Param.r^2; % Cross section area
% Param.J1 = pi*Param.r^4/2; % Polar inertia moment
% Param.J2 = pi*Param.r^4/4; % Inertia moment
% Param.J3 = pi*Param.r^4/4; % Inertia moment
% %Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
% Param.H = diag([1.5708e-2 , 0.7854e-2 , 0.7854e-2 , 3.14e2 ,3.14e2 , 3.14e2]); %Hooke tensor
% 
% Param.E = 2.638984243265842e+09; % Young modulu
% Param.E = 3.000205601e+09;
Param.E = 3e+09;
Param.G=8.543*10^4; % shear modulu
%Param.r = 0.02; % rod radius
Param.rho = 1240*0.20; % Mass density
%mp=0;
dis=0.043489875000000;
%u=2.788212272667389e-04;
%u=0.001116146403359
%u=2.478604539739636e-04;
%u= 1e-05;
%u=2.42e-4;
%u=2.14e-04;
u=1.696782934742007e-04;
%u = 1.5e-4;
%u=  1.327234111222623e-04;

%Platform parameters
Param.a=0.014481;
Param.c=0.08;
Param.cl=0.03;
%Param.b=0.0402;
Param.b=dis;
Param.rhop=1240;
Param.Vp=Param.a*Param.b*Param.c;
Param.mp=Param.rhop*Param.Vp;%Param.Vp*Param.rhop;%Param.Vp*Param.rhop;%Param.Vp*Param.rhop;
Param.Jp1=1/12*Param.mp*(Param.b^2+Param.c^2);
Param.Jp2=1/12*Param.mp*(Param.a^2+Param.c^2);
Param.Jp3=1/12*Param.mp*(Param.a^2+Param.b^2);
Param.Mp =diag([Param.Jp1,Param.Jp2,Param.Jp3,Param.mp,Param.mp,Param.mp]);%cross-sectional inertia


Param.g =[0 ; 0 ; 0 ; 9.81 ; 0; 0]; % gravity
Param.L = 0.20;%4121500000000; % Rod length
Param.d=0.001;
Param.A = Param.cl*Param.d; % Cross section area
%Param.rho*Param.A
Param.J1 = Param.cl*Param.d*(Param.cl^2+Param.d^2)/12; % Polar inertia moment
Param.J2 = Param.cl^3*Param.d/12; % Inertia moment
Param.J3 = Param.cl*Param.d^3/12; % Inertia moment
%Param.mu=1000;
%Param.E*Param.J3/30.90

%Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A, Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.J3/u, Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia
Param.dX=Param.L/100; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);



Param.Fep=[0;Param.mp*9.81;0];
Param.Fep=0*[0;40;0];
Param.mu=10^-3;
%Tendons parameters
Param.Rb= 0.016;
Param.c1=pi/2;
Param.c2=-pi/2;
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)]; % coordinate of the tendon in a cross-section
Param.n=3;  %strain modes number
Param.na=2; %number of actuated strains
%Param.dX=Param.L/500;%length of segments
Param.n_lambda=6; %dimension of Lambda
Y=0:Param.dX:Param.L;
Param.Tendons_list=zeros(3*(Param.n_X+1),2);
for i=1:Param.n_X+1
%Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi,Param.Rb)]; % Parallel case
%Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-i*Param.dX))]; % Convergent routing
%Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-2*i*Param.dX))]; % Convergent routing
 %if (i-1)*Param.dX>Param.L/1000
     Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi,Param.Rb)];
 %end
end
%T=0;
%for i=1:100
%T=current2tension(2000)
T1=0
T2=0
Param.Forces_Tendons=[0
T1
T2
0]; 

%  Param.E = 3*10^9; % Young modulu
% Param.G=8.543*10^4; % shear modulu
% %Param.r = 0.02; % rod radius
% Param.rho = 1240*0.2; % Mass density
% %mp=0;
% dis=0.1;
% %u=2.788212272667389e-04;
% %u=0.001116146403359
% %u=2.478604539739636e-04;
% %u= 1.5e-4;
% u=1.696782934742007e-04;
% %Platform parameters
% Param.a=0.014481;
% Param.c=0.08;
% Param.cl=0.03;
% Param.b=dis/5;%0.0402;
% Param.rhop=1240;
% Param.Vp=Param.a*Param.b*Param.c;
% Param.mp=Param.Vp*Param.rhop;%Param.Vp*Param.rhop;%Param.Vp*Param.rhop;
% Param.Jp1=1/12*Param.mp*(Param.b^2+Param.c^2);
% Param.Jp2=1/12*Param.mp*(Param.a^2+Param.c^2);
% Param.Jp3=1/12*Param.mp*(Param.a^2+Param.b^2);
% %Param.Mp =diag([Param.Jp1,Param.Jp2,Param.Jp3,Param.mp,Param.mp,Param.mp]);%cross-sectional inertia
% Param.Mp =diag([Param.Jp3,Param.mp,Param.mp]);%cross-sectional inertia
% 
% Param.g =[0 ; 0 ; 0 ; 9.81 ; 0; 0]; % gravity
% Param.L = 0.20;%4121500000000; % Rod length
% Param.d=0.001;
% Param.A = Param.cl*Param.d; % Cross section area
% Param.J1 = Param.cl*Param.d*(Param.cl^2+Param.d^2)/12; % Polar inertia moment
% Param.J2 = Param.cl^3*Param.d/12; % Inertia moment
% Param.J3 = Param.cl*Param.d^3/12; % Inertia moment
% %Param.mu=1000;
% %Param.J3/Param.A
% %Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A, Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
% Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.J3/u, Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
% Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia
% Param.dX=Param.L/100; %spatial step
% Param.n_X=length(Param.dX:Param.dX:Param.L);
% 
% %Param.mu=0.4;
% %Param.D = Param.mu*diag([Param.J1 , 3*Param.J2 , 3*Param.J3 , 3*Param.A , Param.A , Param.A]);%Damping matrix
% Param.D=10^(-5)*eye(6);%Damping matrix
% Param.duree=1;%simulation time
% Param.dt=0.001; %time step
% Param.dX=Param.L/100; %spatial step
% Param.n_X=length(Param.dX:Param.dX:Param.L);
% Param.n_t=length(Param.dt:Param.dt:Param.duree);
% Param.n_lambda=6;
% 
% 
% 
% %Tendons parameters
% Param.Rb= 0.016; % Distance between a tendon and the backbone
% Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)];% coordinate of the tendon in a cross-section
% %Param.Tendons_list = [Param.Tendon_coordinate(0),Param.Tendon_coordinate(pi/2),Param.Tendon_coordinate(pi),Param.Tendon_coordinate(3*pi/2)]; % 4 tendons separated with an angle of pi/2
% Param.n=3;  %strain modes number
% Param.na=2; %number of actuated strains
% 
% Y=0:Param.dX:Param.L;
% Param.Tendons_list=zeros(3*(Param.n_X+1),2);
% for i=1:Param.n_X+1
% %Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi,Param.Rb)]; % Parallel case
% %Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-i*Param.dX))]; % Convergent routing
% %Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-2*i*Param.dX))]; % Convergent routing
% %Parallel truncated
% %if (i-1)*Param.dX>Param.L/50
% Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi,Param.Rb)];
% %end
% end
% %Param.Forces_Tendons = zeros(4,1); % tension in each tendon
% 
% %T=current2tension(1000) 
%  T1= 15.746446056012470;
%  T2= 24.903293743502662;
% Param.Forces_Tendons=[0
% T1
% T2
% 0];   
Param.Ftip=0*[0;-1/sqrt(2);1/sqrt(2);0;0;0];
%Param.Forces_Tendons(1)=1;

%Initial conditions
Param.theta01=-0*pi/70;
Param.theta02=-0*pi/70;
R01=[cos(Param.theta01) -sin(Param.theta01) 0; sin(Param.theta01) cos(Param.theta01) 0; 0 0 1];
R02=[cos(Param.theta02) -sin(Param.theta02) 0; sin(Param.theta02) cos(Param.theta02) 0; 0 0 1];
Param.g1=eye(4);
Param.g1(1:3,1:3)=R01;
Param.y1=0;
Param.g1(2,4)=Param.y1;
Param.g2=eye(4);
Param.y2=Param.y1+dis;%0.027175500000000;
%Param.y2=0.03;
Param.g2(1:3,1:3)=R02;
Param.g2(2,4)=Param.y2;
Param.CM1=[Param.a/2;(Param.y2-Param.y1)/2];
Param.CM2=[Param.a/2;-(Param.y2-Param.y1)/2];

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

% Param.B     = [0 0 0 ;0 0 0;1 0 0;0 1 0;0 0 1;0 0 0];
% Param.B_bar = [1 0 0;0 1 0;0 0 0;0 0 0;0 0 0;0 0 1];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;0];

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
Param.F_p=[0;0;0];
 init(14:15)=[Param.L+Param.a/2;(Param.y2+Param.y1)/2];
 %   init=[    5.767963539474240
 %   4.255950051640640
 %   7.814798085483162
 %   0.000054135438839
 %  -0.000004252186027
 %  -0.000005885238970
 %   5.767963539473167
 %   5.448184480654892
 % -31.423718309334678
 %  -0.000175639181651
 %   0.000002192331196
 %   0.000009510832515
 %   1.177365369622749
 %   0.161013203013013
 %   0.082905876680758
 %   0.079884624869980
 %   4.714005216032118
 %   1.732137215995682
 %  -0.153423344286255
 %  -4.071101692713156
 %  -1.730940563502999];

%    init=[      4.684127300946545
%    3.062253503415598
%    6.936878367313686
%    0.000061714775343
%   -0.000002788143934
%   -0.000004151590065
%    4.684127300945152
%    3.741300810535179
%  -30.047016201023222
%   -0.000192366485766
%    0.000000912903330
%    0.000007683629476
%    0.956131090860054
%    0.172519122467663
%    0.074885910922101
%    0.065923414411091
%    5.365503475379665
%    1.745646894675108
%   -0.161909584478044
%   -4.668612532970146
%   -1.745646894675108
% ];
 % init=[    -8.036678421183611
 %  -0.093494162987599
 %  32.548365065510566
 %  -0.009079891201711
 %  -0.000007519328411
 %   0.001035596544133
 %  -8.036769033264822
 %  -0.066124669113998
 % -10.899219074118406
 %   0.003081420009149
 %  -0.000006499143455
 %  -0.000569841669149
 %  -1.640462641086406
 %   0.109234999808697
 %  -0.119107798384230
 %   0.150632775056447
 %  -3.361112951340991
 %   3.589711378265617
 %  -0.077969983902259
 %   3.361112917481296
 %  -3.589711368027869];
%for t=[0:0.001:10]
  %    init=[  -1.504989513688684 %T=4N
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

  %    init=[  -3.900118563682340 % T=3.45N
  % -0.000000000000013
  % 24.702340514112308
  % -0.000109207830704
  %  0.000000000000000
  %  0.000005590326500
  % -3.900118563682343
  % -0.000000000000010
  % -5.550028565767714
  %  0.000059494047545
  %  0.000000000000000
  % -0.000002418871729
  % -0.780023712736468
  %  0.177599049404503
  % -0.060116111413517
  %  0.083293388749573
  % -3.347675268091437
  %  1.376124771260256
  % -0.026382909192019
  %  3.347675268091437
  % -1.376124771260256];

  %    init=[  -5.412893507897341  % T=3.45N
  %  0.000000254952119
  % 26.572746816567133
  % -0.000102395883676
  % -0.000000000000011
  %  0.000008805144415
  % -5.412893507897357
  %  0.000000213864259
  % -7.357494759017236
  %  0.000054741356679
  %  0.000000000000171
  % -0.000004359336601
  % -1.082578701579470
  %  0.156452304972675
  % -0.083068400030528
  %  0.084390725197994
  % -2.893353828278781
  %  1.739436938433977
  % -0.035203709866058
  %  2.893353814890353
  % -1.739436950413657];
  
 %     init=[ -10.104343199150291%T=4N
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
 %   0.076117184007060
 %  -0.120248612677159
 %   0.083392714504166
 %  -1.478386599850913
 %   2.356167496296465
 %  -0.058260142305702
 %   1.478386599838229
 %  -2.356167496297970];


%       for T=[0:0.1:20]
%      %     T
% % % % %    Param.g =[0 ; 0 ; 0 ; -g ; 0; 0]; % gravity
%         Param.Forces_Tendons=[0
%                               T
%                               T
%                              0]; 
%Param.F_p=Param.F_p+[0;1;0];

 %    init=[   -5.235987756315321
 %  15.213413762917051
 %   0.329374441147704
 %  -0.000214813185683
 %   0.020595413547377
 %  -0.019216299349561
 %   5.235987756315319
 % -15.213413762917092
 %  -0.329374441147730
 %  -0.000214813185684
 %   0.020595413547377
 %  -0.019216299349561
 %  -0.000000000000000
 %   0.196238916869977
 %   0.050000000000000
 %   0.025469774413990
 %   0.000000000000003
 %  -0.792797720563080
 %  -0.025469774413990
 %  -0.000000000000003
 %   0.792797720563080];

 options=optimset('fsolve');
 %options.MaxFunEvals =1;
 options.MaxIter = 150;
 options.TolFun=1e-20;
 %options.StepTolerance=10^-5;
 %options.exitflag=4;
%options.OptimalityTolerance=10^(-5);
 options.Display = 'iter';
options.Jacobian='on';
% options.Algorithm='levenberg-marquardt';
 %options.Algorithm='trust-region-dogleg';
tic
%  [x,f,a,b,J]=fsolve(@(var)Static(var,Param),init,options);
% x
% lb=zeros(2*Param.n*Param.na+Param.n_lambda+3,1);
% ub=zeros(2*Param.n*Param.na+Param.n_lambda+3,1);
% lb(1:2*Param.n*Param.na)=-10;
% lb(2*Param.n*Param.na+1)=-pi;
% lb(2*Param.n*Param.na+2)=0;
% lb(2*Param.n*Param.na+3)=-0.2;
% lb(2*Param.n*Param.na+4:2*Param.n*Param.na+Param.n_lambda+3)=-inf;
% ub(1:2*Param.n*Param.na)=10;
% ub(2*Param.n*Param.na+1)=pi;
% ub(2*Param.n*Param.na+2)=0.22;
% ub(2*Param.n*Param.na+3)=0.2;
% ub(2*Param.n*Param.na+4:2*Param.n*Param.na+Param.n_lambda+3)=inf;
% A            =  []; 
% b            =  [];
% Aeq          =  [];
% beq          =  [];
% C            =  [];
% options  =  optimoptions( 'fmincon', ...
%                           'ConstraintTolerance', 1e-7, ...
%                           'FiniteDifferenceStepSize', 1e-7, ...
%                           'OptimalityTolerance', 1e-7, ...
%                           'StepTolerance', 1e-7, ...
%                           'FunctionTolerance', 1e-4, ...
%                           'SpecifyObjectiveGradient', false, ...
%                           'Display', 'off');

% options  =  optimoptions( 'fsolve', ...
%                           'FiniteDifferenceStepSize', 1e-4, ...
%                           'OptimalityTolerance', 1e-4, ...
%                           'StepTolerance', 1e-4, ...
%                           'FunctionTolerance', 1e-4, ...
%                           'SpecifyObjectiveGradient', false, ...
%                           'Display', 'off');
[x,f,a,b,J]=fsolve(@(var)Static_gradient(var,Param),init,options);
%x
%[x,f,a,b,J]=fsolve(@(var)Static_pivot(var,Param),init,options);
%x=fmincon(@(var)Static(var,Param),init,A,b,Aeq,beq,lb,ub,C,options);


%options.Jacobian='off';
%[x,f,a,b,J1]=fsolve(@(var)Static(var,Param),init,options);
%[x,f,a,b,J2]=fsolve(@(var)Static_V2(var,Param),init,options);
%J
%J1
%J-J1
init=x;
q1=x(1:Param.n*Param.na);
q2=x(Param.n*Param.na+1:2*Param.n*Param.na); 
qp=x(2*Param.n*Param.na+1:2*Param.n*Param.na+3);
y(length(x)*(i-1)+1:length(x)*i)=x;
%norm(f)
toc

r1=zeros(3,Param.n_X+1);
r1(2,1)=Param.y1;
r2=zeros(3,Param.n_X+1);
r2(2,1)=Param.y2;
rm1=zeros(3,Param.n_X+1);
rm1(2,1)=Param.y1-Param.Rb;
rm2=zeros(3,Param.n_X+1);
rm2(2,1)=Param.y2+Param.Rb;
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
    rm1(:,j+1)=r1(:,j+1)+R1*Param.Tendon_coordinate(pi, Param.Rb);
    rm2(:,j+1)=r2(:,j+1)+R2*Param.Tendon_coordinate(0, Param.Rb);
end  
g1=[R1 r1(:,end);0 0 0 1];
g2=[R2 r2(:,end);0 0 0 1];
Rp=[cos(qp(1)) -sin(qp(1)); sin(qp(1)) cos(qp(1))];

    figure
    hold on
    l_rec=0.02;
    plot(r1(1,:),r1(2,:),'b','LineWidth', 4); %title('Planar TACR'); axis([-0.1 0.25 -0.1 0.1 ]);grid on; daspect([1 1 1])
    axis([-Param.L 2*Param.L -Param.L Param.L ]);grid on; daspect([1 1 1])
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow;
    %hold on    
    plot(rm1(1,:),rm1(2,:),'g','LineWidth', 2);
    plot(rm2(1,:),rm2(2,:),'g','LineWidth', 2);
    plot(r2(1,:),r2(2,:),'b','LineWidth', 4);drawnow;
    
   
   plot([-0.1;rm1(1,1)],[rm1(2,1);rm1(2,1)],'g','LineWidth', 2);
   plot([-0.1;rm2(1,1)],[rm2(2,1);rm2(2,1)],'g','LineWidth', 2);
    
     plot([r1(1,1);rm1(1,1)],[r1(2,1);rm1(2,1)],'b','LineWidth', 2);drawnow;
     plot([r2(1,1);rm2(1,1)],[r2(2,1);rm2(2,1)],'b','LineWidth', 2);drawnow;

    plot([r1(1,20);rm1(1,20)],[r1(2,20);rm1(2,20)],'b','LineWidth', 2);drawnow;
    plot([r2(1,20);rm2(1,20)],[r2(2,20);rm2(2,20)],'b','LineWidth', 2);drawnow;

    plot([r1(1,40);rm1(1,40)],[r1(2,40);rm1(2,40)],'b','LineWidth', 2);drawnow;
    plot([r2(1,40);rm2(1,40)],[r2(2,40);rm2(2,40)],'b','LineWidth', 2);drawnow;

    plot([r1(1,60);rm1(1,60)],[r1(2,60);rm1(2,60)],'b','LineWidth', 2);drawnow;
    plot([r2(1,60);rm2(1,60)],[r2(2,60);rm2(2,60)],'b','LineWidth', 2);drawnow;

    plot([r1(1,80);rm1(1,80)],[r1(2,80);rm1(2,80)],'b','LineWidth', 2);drawnow;
    plot([r2(1,80);rm2(1,80)],[r2(2,80);rm2(2,80)],'b','LineWidth', 2);drawnow;
    % 
     plot([r1(1,end);rm1(1,end)],[r1(2,end);rm1(2,end)],'b','LineWidth', 2);drawnow;
     plot([r2(1,end);rm2(1,end)],[r2(2,end);rm2(2,end)],'b','LineWidth', 2);drawnow;

    plot(qp(2),qp(3),'r*','LineWidth', 4);drawnow;
%       xrec=[qp(2); qp(2)+2*Param.a/2*cos(qp(1))];
%       yrec=[qp(3); qp(3)+2*Param.a/2*sin(qp(1))];
%       plot(xrec,yrec,'r','LineWidth', 4);
     xrec=[qp(2)+(Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1))); qp(2)+(Param.a/2*cos(qp(1))+Param.b/2*sin(qp(1))) ;qp(2)-(Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1)));qp(2)+(-Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1)));qp(2)+(Param.a/2*cos(qp(1))-Param.b/2*sin(qp(1)))];
     yrec=[qp(3)+(+Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1))) ;qp(3)+(Param.a/2*sin(qp(1))-Param.b/2*cos(qp(1))) ;qp(3)-(Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1)));qp(3)+(-Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1)));qp(3)+(Param.a/2*sin(qp(1))+Param.b/2*cos(qp(1)))];
     plot(xrec,yrec,'r','LineWidth', 4);drawnow;
     hold off
%end
 %end
%end
   % quiver3(r(1,end),r(2,end),r(3,end),r(1,end)+Ft(1),r(2,end)+Ft(2),r(3,end)+Ft(3));