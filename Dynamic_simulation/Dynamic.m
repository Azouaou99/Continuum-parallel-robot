function res = Dynamic(var,q0,qd0,i,Param)

q1=var(1:Param.n*Param.na);
q2=var(Param.n*Param.na+1:2*Param.n*Param.na);
qp=var(2*Param.n*Param.na+1:2*Param.n*Param.na+3);
qd1=(q1-q0(1:Param.n*Param.na))/Param.dt;%var(2*Param.n*Param.na+1:3*Param.n*Param.na);
qd2=(q2-q0(Param.n*Param.na+1:2*Param.n*Param.na))/Param.dt; %var(3*Param.n*Param.na+1:4*Param.n*Param.na);
qdp=(qp-q0(2*Param.n*Param.na+1:2*Param.n*Param.na+3))/Param.dt;
qdd1=(qd1-qd0(1:Param.n*Param.na))/Param.dt;%var(2*Param.n*Param.na+1:3*Param.n*Param.na);
qdd2=(qd2-qd0(Param.n*Param.na+1:2*Param.n*Param.na))/Param.dt; %var(3*Param.n*Param.na+1:4*Param.n*Param.na);
qddp=(qdp-qd0(2*Param.n*Param.na+1:2*Param.n*Param.na+3))/Param.dt;
lambda=var(2*Param.n*Param.na+4:2*Param.n*Param.na+Param.n_lambda+3);
Rp=[cos(qp(1)) -sin(qp(1)) 0; sin(qp(1)) cos(qp(1)) 0; 0 0 1];
ql=[q1;q2];
SM1=zeros(Param.n*Param.na,Param.n*Param.na);
SM2=zeros(Param.n*Param.na,Param.n*Param.na);
SC1=zeros(Param.n*Param.na,Param.n*Param.na);
SC2=zeros(Param.n*Param.na,Param.n*Param.na);
SFe1=zeros(Param.n*Param.na,1);
SFe2=zeros(Param.n*Param.na,1);
SFa1=zeros(Param.n*Param.na,1);
SFa2=zeros(Param.n*Param.na,1);

g1=Param.g1;
g2=Param.g2;
J1= zeros(6,Param.n*Param.na);
J2= zeros(6,Param.n*Param.na);
I1= zeros(6,Param.n*Param.na);
I2= zeros(6,Param.n*Param.na);

Jd1= zeros(6,Param.n*Param.na);
Jd2= zeros(6,Param.n*Param.na);
Id1= zeros(6,Param.n*Param.na);
Id2= zeros(6,Param.n*Param.na);
% Y1=0:Param.dX:Param.L;
%  for j=1:Param.n_X
%      Jp0=Jp(:,(j-1)*Param.n*Param.na+1:(j)*Param.n*Param.na);
%      Jp(:,(j)*Param.n*Param.na+1:(j+1)*Param.n*Param.na)=Jacobian_subseg(Jp0,Y1(j+1),q,Param);
%      Jp1=Jp(:,(j)*Param.n*Param.na+1:(j+1)*Param.n*Param.na);
%      Jdp0=Jdp(:,(j-1)*Param.n*Param.na+1:(j)*Param.n*Param.na);
%      Jdp(:,(j)*Param.n*Param.na+1:(j+1)*Param.n*Param.na)=dJacobian_subseg(Jp1,Jp0,Jdp0,Y1(j+1),q,qd,Param);
%  end
Y=0:Param.dX:Param.L;
phival=Phi(Param.na,Param.n,Y(1),Param.L);
H1=Ftendons(Y(1),1,q1,Param);
H2=Ftendons(Y(1),1,q2,Param);
fFa10=-phival'*Param.B'*H1;
fFa20=-phival'*Param.B'*H2;
fM10=J1'*Param.M*J1;
fM20=J2'*Param.M*J2;
fC10=J1'*(Param.M*Jd1 -ad_func(J1*qd1)' * Param.M*J1);
fC20=J2'*(Param.M*Jd2 -ad_func(J2*qd2)' * Param.M*J2);
fFe10=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
fFe20=-J2'*Param.M*Ad_function(Inverse_T(g2))*Param.g;
P1=Ad_function(g1)*Param.B*phival;
P2=Ad_function(g2)*Param.B*phival;
Pd1=Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1;
Pd2=Ad_function(g2)*ad_func(Param.B*phival*qd2)*J2;
h=Param.dX;
for k=2:Param.n_X+1
%     xiaval1=Param.xi_a0+phival*q1;
%     xival1=Param.B*xiaval1+Param.B_bar*Param.xi_c;
%     xiaval2=Param.xi_a0+phival*q2;
%     xival2=Param.B*xiaval2+Param.B_bar*Param.xi_c;
%     phival=Phi(Param.na,Param.n,Y(k),Param.L);
%     g1=g1*expm(Hat(xival1)*(Param.dX));
%     g2=g2*expm(Hat(xival2)*(Param.dX));
     %g1=g1*expm2(Hat(xival1)*(Param.dX),xival1*(Param.dX));
     %g2=g2*expm2(Hat(xival2)*(Param.dX),xival2*(Param.dX));
    x1=Y(k-1)+h/2-sqrt(3)*h/6;
    x2=Y(k-1)+h/2+sqrt(3)*h/6;
    phival1=Phi(Param.na,Param.n,x1,Param.L);
    phival2=Phi(Param.na,Param.n,x2,Param.L);
    xiaval11=Param.xi_a0+phival1*q1;
    xival11=Param.B*xiaval11+Param.B_bar*Param.xi_c;
    xiaval21=Param.xi_a0+phival1*q2;
    xival21=Param.B*xiaval21+Param.B_bar*Param.xi_c;
   
    xiaval12=Param.xi_a0+phival2*q1;
    xival12=Param.B*xiaval12+Param.B_bar*Param.xi_c;
    xiaval22=Param.xi_a0+phival2*q2;
    xival22=Param.B*xiaval22+Param.B_bar*Param.xi_c;
    
    om_h1=h/2*(xival11+xival12)+sqrt(3)*h^2/12* ad_func(xival11)*xival12;
    om_h2=h/2*(xival21+xival22)+sqrt(3)*h^2/12* ad_func(xival21)*xival22;    
    g1=g1*expm(Hat(om_h1));
    g2=g2*expm(Hat(om_h2));
    phival=Phi(Param.na,Param.n,Y(k),Param.L);
    
    N1=Ad_function(g1)*Param.B*phival;
    N2=Ad_function(g2)*Param.B*phival;
    I1=I1+(P1+N1)*Param.dX/2;
    I2=I2+(P2+N2)*Param.dX/2;
    J1=Ad_function(Inverse_T(g1))*I1;
    J2=Ad_function(Inverse_T(g2))*I2;
    Nd1=Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1;
    Nd2=Ad_function(g2)*ad_func(Param.B*phival*qd2)*J2;
    Id1=Id1+(Pd1+Nd1)*Param.dX/2;
    Id2=Id2+(Pd2+Nd2)*Param.dX/2;
    Jd1=-Ad_function(Inverse_T(g1))*Id1;
    Jd2=-Ad_function(Inverse_T(g2))*Id2;
    
    H1=Ftendons(Y(k),k,q1,Param);
    H2=Ftendons(Y(k),k,q2,Param);
    fFa1=-phival'*Param.B'*H1;
    fFa2=-phival'*Param.B'*H2;
    fM1=J1'*Param.M*J1;
    fM2=J2'*Param.M*J2;
    fC1=J1'*(Param.M*Jd1 -ad_func(J1*qd1)' * Param.M*J1);
    fC2=J2'*(Param.M*Jd2 -ad_func(J2*qd2)' * Param.M*J2);
    fFe1=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
    fFe2=-J2'*Param.M*Ad_function(Inverse_T(g2))*Param.g;
    %fFa=-phival'*Param.B'*(Ftendonsval);
 
    if k==Param.n_X+1
        fM1n=fM1;
        fM2n=fM2;
        fC1n=fC1;
        fC2n=fC2;
        fFe1n=fFe1;
        fFe2n=fFe2;
        fFa1n=fFa1;
        fFa2n=fFa2;
    end
    if and(k~=Param.n_X+1,k~=1)
        SM1=SM1+fM1;
        SM2=SM2+fM2;
        SC1=SC1+fC1;
        SC2=SC2+fC2;
        SFe1=SFe1+fFe1;
        SFe2=SFe2+fFe2;
        SFa1=SFa1+fFa1;
        SFa2=SFa2+fFa2;
    end
    %Q1=Q1 + Q1_X*Param.dX;
    %r1=r1 + r1_X*Param.dX;
P1=N1;
P2=N2;
Pd1=Nd1;
Pd2=Nd2;
end
r1=g1(1:3,4);
r2=g2(1:3,4);
R1=g1(1:3,1:3);
R2=g2(1:3,1:3);

Mq1=Param.L/(Param.n_X)*((fM10+fM1n)/2 + SM1);
Mq2=Param.L/(Param.n_X)*((fM20+fM2n)/2 + SM2);
Mq=[Mq1, zeros(Param.n*Param.na,3+Param.n*Param.na);zeros(Param.n*Param.na), Mq2, zeros(Param.n*Param.na,3);zeros(3,2*Param.n*Param.na) Param.Mp];
Cqqd1=Param.L/(Param.n_X)*((fC10+fC1n)/2 + SC1);
Cqqd2=Param.L/(Param.n_X)*((fC20+fC2n)/2 + SC2);
Cqqd=[Cqqd1, zeros(Param.n*Param.na,3+Param.n*Param.na);zeros(Param.n*Param.na), Cqqd2 , zeros(Param.n*Param.na,3); zeros(3,2*Param.n*Param.na+3)];
Fe1=Param.L/(Param.n_X)*((fFe10+fFe1n)/2 + SFe1);
Fe2=Param.L/(Param.n_X)*((fFe20+fFe2n)/2 + SFe2);
F_p=[0;-Param.mp*Param.g(4:5)];
Fe=[Fe1;Fe2;F_p];
H1q= Param.dX.*((fFa10+fFa1n)./2 + SFa1);
H2q= Param.dX.*((fFa20+fFa2n)./2 + SFa2);
Hq=[H1q, zeros(Param.n*Param.na,2);zeros(Param.n*Param.na,2), H2q; zeros(3,4)];
Param.Keps;
Keps=[Param.Keps, zeros(Param.n*Param.na,Param.n*Param.na+3);zeros(Param.n*Param.na), Param.Keps zeros(Param.n*Param.na,3);zeros(3,2*Param.n*Param.na+3)];
Deps=[Param.Deps, zeros(Param.n*Param.na,Param.n*Param.na+3);zeros(Param.n*Param.na), Param.Deps zeros(Param.n*Param.na,3);zeros(3,2*Param.n*Param.na+3)];
Rc1=[cos(Param.c1) -sin(Param.c1) 0; sin(Param.c1) cos(Param.c1) 0; 0 0 1];
Rc2=[cos(Param.c2) -sin(Param.c2) 0; sin(Param.c2) cos(Param.c2) 0; 0 0 1];
Psi1=Hat_inv(Rp*Rc1*R1'- R1*Rc1'*Rp');
Psi1=Psi1(3);
Psi2=qp(2:3)-r1(1:2)-Rp(1:2,1:2)*Param.CM1;
Psi3=Hat_inv(Rp*Rc2*R2'- R2*Rc2'*Rp');
Psi3=Psi3(3);
Psi4=qp(2:3)-r2(1:2)-Rp(1:2,1:2)*Param.CM2;
Psi=[Psi1;Psi2;Psi3;Psi4];
P_p=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_p=[0;0;0;1];
P_R=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_R=[1 0 0;0 1 0;0 0 1;0 0 0];

Grad_g1_q=zeros(4,2*4*Param.n*Param.na);
Grad_g2_q=zeros(4,2*4*Param.n*Param.na);
Psi1_q=zeros(3,2*Param.n*Param.na);
Psi2_q=zeros(3,2*Param.n*Param.na);
Psi3_q=zeros(3,2*Param.n*Param.na);
Psi4_q=zeros(3,2*Param.n*Param.na);
for i=1:Param.n*Param.na
Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J1(:,i));
Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))=g2*Hat(J2(:,i));
Psi1_q(:,i)=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'-(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)*Rc1'*Rp');
Psi2_q(:,i)=-P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p);
Psi3_q(:,Param.n*Param.na+i)=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)'- (P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)*Rc2'*Rp');
Psi4_q(:,Param.n*Param.na+i)=-P_p*(Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_p);
end
Rp_qp1=[-sin(qp(1)) -cos(qp(1)) 0; cos(qp(1)) -sin(qp(1)) 0; 0 0 0];
Psi1_q=Psi1_q(3,:);
Psi2_q=Psi2_q(1:2,:);
Psi3_q=Psi3_q(3,:);
Psi4_q=Psi4_q(1:2,:);
Psi1_qp1=Hat_inv(Rp_qp1*Rc1*(P_R*g1*Q_R)'- (P_R*g1*Q_R)*Rc1'*Rp_qp1');
Psi1_qp1=Psi1_qp1(3);
Psi1_qp2=0;
Psi1_qp3=0;

Psi2_qp1=-Rp_qp1(1:2,1:2)*Param.CM1;
Psi2_qp2=[1;0];
Psi2_qp3=[0;1];

Psi3_qp1=Hat_inv(Rp_qp1*Rc2*(P_R*g2*Q_R)'- (P_R*g2*Q_R)*Rc2'*Rp_qp1');
Psi3_qp1=Psi3_qp1(3);
Psi3_qp2=0;
Psi3_qp3=0;

Psi4_qp1=-Rp_qp1(1:2,1:2)*Param.CM2;
Psi4_qp2=[1;0];
Psi4_qp3=[0;1];


% Psi_q=[Psi1_q1 Psi1_q2 Psi1_q3 Psi1_q4 Psi1_q5 Psi1_q6 Psi1_q7 Psi1_q8 Psi1_q9 Psi1_q10 Psi1_q11 Psi1_q12 Psi1_qp1 Psi1_qp2 Psi1_qp3
%     Psi2_q1 Psi2_q2 Psi2_q3 Psi2_q4 Psi2_q5 Psi2_q6 Psi2_q7 Psi2_q8 Psi2_q9 Psi2_q10 Psi2_q11 Psi2_q12 Psi2_qp1 Psi2_qp2 Psi2_qp3
%     Psi3_q1 Psi3_q2 Psi3_q3 Psi3_q4 Psi3_q5 Psi3_q6 Psi3_q7 Psi3_q8 Psi3_q9 Psi3_q10 Psi3_q11 Psi3_q12 Psi3_qp1 Psi3_qp2 Psi3_qp3
%     Psi4_q1 Psi4_q2 Psi4_q3 Psi4_q4 Psi4_q5 Psi4_q6 Psi4_q7 Psi4_q8 Psi4_q9 Psi4_q10 Psi4_q11 Psi4_q12 Psi4_qp1 Psi4_qp2 Psi4_qp3];
Psi_q=[Psi1_q Psi1_qp1 Psi1_qp2 Psi1_qp3
    Psi2_q Psi2_qp1 Psi2_qp2 Psi2_qp3
    Psi3_q Psi3_qp1 Psi3_qp2 Psi3_qp3
    Psi4_q Psi4_qp1 Psi4_qp2 Psi4_qp3];
qdd=[qdd1;qdd2;qddp];
qd=[qd1;qd2;qdp];
q=[q1;q2;qp];
res1=Mq*qdd-Hq*Param.Forces_Tendons+Fe+Cqqd*qd+Keps*q+Deps*qd+Psi_q'*lambda;%2nd derivative of q
%res1=-Fa+Keps*q+Deps*qd+Psi_q'*lambda;%2nd derivative of q
res2=Psi;
res=[res1;res2]; %Objective function

%Wt= ql'*[H1q, zeros(Param.n*Param.na,2);zeros(Param.n*Param.na,2), H2q]*Param.Forces_Tendons
%El=ql'*[Param.Keps, zeros(Param.n*Param.na,Param.n*Param.na);zeros(Param.n*Param.na), Param.Keps]*ql
%E=ql'*[Param.Keps, zeros(Param.n*Param.na,Param.n*Param.na);zeros(Param.n*Param.na), Param.Keps]*ql + ql'*[H1q, zeros(Param.n*Param.na,2);zeros(Param.n*Param.na,2), H2q]*Param.Forces_Tendons
end