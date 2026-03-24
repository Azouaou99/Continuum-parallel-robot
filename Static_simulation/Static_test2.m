function res = Static_test2(var,Param)
q1=var(1:Param.n*Param.na);
q2=var(Param.n*Param.na+1:2*Param.n*Param.na);
qp=var(2*Param.n*Param.na+1:2*Param.n*Param.na+3);
Rp=[cos(qp(1)) -sin(qp(1)) 0; sin(qp(1)) cos(qp(1)) 0; 0 0 1];
rp=[qp(2);qp(3);0];
lambda=var(2*Param.n*Param.na+4 : 2*Param.n*Param.na+3+Param.n_lambda);
q=[q1;q2;qp];
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
Y=0:Param.dX:Param.L;
for k=1:Param.n_X+1
    phival=Phi(Param.na,Param.n,Y(k),Param.L);
    xiaval1=Param.xi_a0+phival*q1;
    xival1=Param.B*xiaval1+Param.B_bar*Param.xi_c;
    xiaval2=Param.xi_a0+phival*q2;
    xival2=Param.B*xiaval2+Param.B_bar*Param.xi_c;

    H1=Ftendons(Y(k),k,q1,Param);
    H2=Ftendons(Y(k),k,q2,Param);
    fFa1=-phival'*Param.B'*H1;
    fFa2=-phival'*Param.B'*H2;
    fFe1=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
    fFe2=-J2'*Param.M*Ad_function(Inverse_T(g2))*Param.g;
    %fFa=-phival'*Param.B'*(Ftendonsval);
    if k==1
        fFe10=fFe1;
        fFe20=fFe2;
        fFa10=fFa1;
        fFa20=fFa2;
    end
    if k==Param.n_X+1
        fFe1n=fFe1;
        fFe2n=fFe2;
        fFa1n=fFa1;
        fFa2n=fFa2;
    end
    if and(k~=Param.n_X+1,k~=1)
        SFe1=SFe1+fFe1;
        SFe2=SFe2+fFe2;
        SFa1=SFa1+fFa1;
        SFa2=SFa2+fFa2;
    end
    %Q1=Q1 + Q1_X*Param.dX;
    %r1=r1 + r1_X*Param.dX;
    P1=Ad_function(g1)*Param.B*phival;
    P2=Ad_function(g2)*Param.B*phival;
    if (k~=Param.n_X+1)
    g1=g1*expm(Hat(xival1)*(Param.dX));
    g2=g2*expm(Hat(xival2)*(Param.dX));
    %g1=g1*expm2(Hat(xival1)*(Param.dX),xival1*(Param.dX));
    %g2=g2*expm2(Hat(xival2)*(Param.dX),xival2*(Param.dX));
    N1=Ad_function(g1)*Param.B*phival;
    N2=Ad_function(g2)*Param.B*phival;
    I1=I1+(P1+N1)*Param.dX/2;
    I2=I2+(P2+N2)*Param.dX/2;
    J1=Ad_function(Inverse_T(g1))*I1;
    J2=Ad_function(Inverse_T(g2))*I2;
    end
end
r1=g1(1:3,4);
r2=g2(1:3,4);
R1=g1(1:3,1:3);
R2=g2(1:3,1:3);

Fe1=Param.L/(Param.n_X)*((fFe10+fFe1n)/2 + SFe1);
Fe2=Param.L/(Param.n_X)*((fFe20+fFe2n)/2 + SFe2);

gp=[Rp,rp;0 0 0 1];
%Fep=Param.Mp*Ad_function(Inverse_T(gp))*Param.g
Fep=-Param.Mp*Param.g;
%Fep=-Param.Mp*Param.g
Fep=[0 0 1 0 0 0 ; 0 0 0 1 0 0;0 0 0 0 1 0]*Fep;
Fe=[Fe1;Fe2;Fep];
H1q= Param.dX.*((fFa10+fFa1n)./2 + SFa1);
H2q= Param.dX.*((fFa20+fFa2n)./2 + SFa2);
Hq=[H1q, zeros(Param.n*Param.na,2);zeros(Param.n*Param.na,2), H2q; zeros(3,4)];
Param.Keps;
Keps=[Param.Keps, zeros(Param.n*Param.na,Param.n*Param.na+3);zeros(Param.n*Param.na), Param.Keps zeros(Param.n*Param.na,3);zeros(3,2*Param.n*Param.na+3)];
Rc1=[cos(Param.c1) -sin(Param.c1) 0; sin(Param.c1) cos(Param.c1) 0; 0 0 1];
Rc2=[cos(Param.c2) -sin(Param.c2) 0; sin(Param.c2) cos(Param.c2) 0; 0 0 1];
%Psi=[X1*(R1'*(r2-r1));X2*(R2'*(r2-r1));sqrt((r1-r2)'*(r1-r2))-Param.y2];
Psi1=Hat_inv(Rp*Rc1*R1'- R1*Rc1'*Rp');
Psi1=Psi1(3);
Psi2=qp(2:3)-r1(1:2)+Rp(1:2,1:2)*Param.CM1;
Psi3=Hat_inv(Rp*Rc2*R2'- R2*Rc2'*Rp');
Psi3=Psi3(3);
Psi4=qp(2:3)-r2(1:2)+Rp(1:2,1:2)*Param.CM2;
Psi=[Psi1;Psi2;Psi3;Psi4];
P_p=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_p=[0;0;0;1];
P_R=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_R=[1 0 0;0 1 0;0 0 1;0 0 0];

Grad_g1_q1=g1*Hat(J1(:,1));
Grad_g1_q2=g1*Hat(J1(:,2));
Grad_g1_q3=g1*Hat(J1(:,3));
Grad_g1_q4=g1*Hat(J1(:,4));
Grad_g1_q5=g1*Hat(J1(:,5));
Grad_g1_q6=g1*Hat(J1(:,6));
Grad_g1_q7=zeros(4);
Grad_g1_q8=zeros(4);
Grad_g1_q9=zeros(4);
Grad_g1_q10=zeros(4);
Grad_g1_q11=zeros(4);
Grad_g1_q12=zeros(4);
%
Grad_g2_q1=zeros(4);
Grad_g2_q2=zeros(4);
Grad_g2_q3=zeros(4);
Grad_g2_q4=zeros(4);
Grad_g2_q5=zeros(4);
Grad_g2_q6=zeros(4);
Grad_g2_q7=g2*Hat(J2(:,1));
Grad_g2_q8=g2*Hat(J2(:,2));
Grad_g2_q9=g2*Hat(J2(:,3));
Grad_g2_q10=g2*Hat(J2(:,4));
Grad_g2_q11=g2*Hat(J2(:,5));
Grad_g2_q12=g2*Hat(J2(:,6));
%[Grad_g1_q1,Grad_g1_q2,Grad_g1_q3,Grad_g1_q4,Grad_g1_q5,Grad_g1_q6,Grad_g2_q7,Grad_g2_q8,Grad_g2_q9,Grad_g2_q10,Grad_g2_q11,Grad_g2_q12]=grad_g(Param.L,q1,q2,Param);
Rp_qp1=[-sin(qp(1)) -cos(qp(1)) 0; cos(qp(1)) -sin(qp(1)) 0; 0 0 0];

Psi1_q1=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q1*Q_R)'-(P_R*Grad_g1_q1*Q_R)*Rc1'*Rp');
Psi1_q1=Psi1_q1(3);
Psi1_q2=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q2*Q_R)'- (P_R*Grad_g1_q2*Q_R)*Rc1'*Rp');
Psi1_q2=Psi1_q2(3);
Psi1_q3=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q3*Q_R)'- (P_R*Grad_g1_q3*Q_R)*Rc1'*Rp');
Psi1_q3=Psi1_q3(3);
Psi1_q4=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q4*Q_R)'-(P_R*Grad_g1_q4*Q_R)*Rc1'*Rp');
Psi1_q4=Psi1_q4(3);
Psi1_q5=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q5*Q_R)'- (P_R*Grad_g1_q5*Q_R)*Rc1'*Rp');
Psi1_q5=Psi1_q5(3);
Psi1_q6=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q6*Q_R)'- (P_R*Grad_g1_q6*Q_R)*Rc1'*Rp');
Psi1_q6=Psi1_q6(3);
Psi1_q7=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q7*Q_R)'- (P_R*Grad_g1_q7*Q_R)*Rc1'*Rp');
Psi1_q7=Psi1_q7(3);
Psi1_q8=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q8*Q_R)'- (P_R*Grad_g1_q8*Q_R)*Rc1'*Rp');
Psi1_q8=Psi1_q8(3);
Psi1_q9=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q9*Q_R)'- (P_R*Grad_g1_q9*Q_R)*Rc1'*Rp');
Psi1_q9=Psi1_q9(3);
Psi1_q10=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q10*Q_R)'- (P_R*Grad_g1_q10*Q_R)*Rc1'*Rp');
Psi1_q10=Psi1_q10(3);
Psi1_q11=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q11*Q_R)'- (P_R*Grad_g1_q11*Q_R)*Rc1'*Rp');
Psi1_q11=Psi1_q11(3);
Psi1_q12=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q12*Q_R)'- (P_R*Grad_g1_q12*Q_R)*Rc1'*Rp');
Psi1_q12=Psi1_q12(3);
Psi1_qp1=Hat_inv(Rp_qp1*Rc1*(P_R*g1*Q_R)'- (P_R*g1*Q_R)*Rc1'*Rp_qp1');
Psi1_qp1=Psi1_qp1(3);
Psi1_qp2=0;
Psi1_qp3=0;

Psi2_q1=-P_p*(Grad_g1_q1*Q_p);
Psi2_q1=Psi2_q1(1:2);
Psi2_q2=-P_p*(Grad_g1_q2*Q_p);
Psi2_q2=Psi2_q2(1:2);
Psi2_q3=-P_p*(Grad_g1_q3*Q_p);
Psi2_q3=Psi2_q3(1:2);
Psi2_q4=-P_p*(Grad_g1_q4*Q_p);
Psi2_q4=Psi2_q4(1:2);
Psi2_q5=-P_p*(Grad_g1_q5*Q_p);
Psi2_q5=Psi2_q5(1:2);
Psi2_q6=-P_p*(Grad_g1_q6*Q_p);
Psi2_q6=Psi2_q6(1:2);
Psi2_q7=-P_p*(Grad_g1_q7*Q_p);
Psi2_q7=Psi2_q7(1:2);
Psi2_q8=-P_p*(Grad_g1_q8*Q_p);
Psi2_q8=Psi2_q8(1:2);
Psi2_q9=-P_p*(Grad_g1_q9*Q_p);
Psi2_q9=Psi2_q9(1:2);
Psi2_q10=-P_p*(Grad_g1_q10*Q_p);
Psi2_q10=Psi2_q10(1:2);
Psi2_q11=-P_p*(Grad_g1_q11*Q_p);
Psi2_q11=Psi2_q11(1:2);
Psi2_q12=-P_p*(Grad_g1_q12*Q_p);
Psi2_q12=Psi2_q12(1:2);
Psi2_qp1=Rp_qp1(1:2,1:2)*Param.CM1;
Psi2_qp2=[1;0];
Psi2_qp3=[0;1];



Psi3_q1=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q1*Q_R)'- (P_R*Grad_g2_q1*Q_R)*Rc2'*Rp');
Psi3_q1=Psi3_q1(3);
Psi3_q2=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q2*Q_R)'- (P_R*Grad_g2_q2*Q_R)*Rc2'*Rp');
Psi3_q2=Psi3_q2(3);
Psi3_q3=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q3*Q_R)'- (P_R*Grad_g2_q3*Q_R)*Rc2'*Rp');
Psi3_q3=Psi3_q3(3);
Psi3_q4=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q4*Q_R)'- (P_R*Grad_g2_q4*Q_R)*Rc2'*Rp');
Psi3_q4=Psi3_q4(3);
Psi3_q5=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q5*Q_R)'- (P_R*Grad_g2_q5*Q_R)*Rc2'*Rp');
Psi3_q5=Psi3_q5(3);
Psi3_q6=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q6*Q_R)'- (P_R*Grad_g2_q6*Q_R)*Rc2'*Rp');
Psi3_q6=Psi3_q6(3);
Psi3_q7=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q7*Q_R)'- (P_R*Grad_g2_q7*Q_R)*Rc2'*Rp');
Psi3_q7=Psi3_q7(3);
Psi3_q8=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q8*Q_R)'- (P_R*Grad_g2_q8*Q_R)*Rc2'*Rp');
Psi3_q8=Psi3_q8(3);
Psi3_q9=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q9*Q_R)'- (P_R*Grad_g2_q9*Q_R)*Rc2'*Rp');
Psi3_q9=Psi3_q9(3);
Psi3_q10=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q10*Q_R)'- (P_R*Grad_g2_q10*Q_R)*Rc2'*Rp');
Psi3_q10=Psi3_q10(3);
Psi3_q11=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q11*Q_R)'- (P_R*Grad_g2_q11*Q_R)*Rc2'*Rp');
Psi3_q11=Psi3_q11(3);
Psi3_q12=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q12*Q_R)'- (P_R*Grad_g2_q12*Q_R)*Rc2'*Rp');
Psi3_q12=Psi3_q12(3);
Psi3_qp1=Hat_inv(Rp_qp1*Rc2*(P_R*g2*Q_R)'- (P_R*g2*Q_R)*Rc2'*Rp_qp1');
Psi3_qp1=Psi3_qp1(3);
Psi3_qp2=0;
Psi3_qp3=0;

Psi4_q1=-P_p*(Grad_g2_q1*Q_p);
Psi4_q1=Psi4_q1(1:2);
Psi4_q2=-P_p*(Grad_g2_q2*Q_p);
Psi4_q2=Psi4_q2(1:2);
Psi4_q3=-P_p*(Grad_g2_q3*Q_p);
Psi4_q3=Psi4_q3(1:2);
Psi4_q4=-P_p*(Grad_g2_q4*Q_p);
Psi4_q4=Psi4_q4(1:2);
Psi4_q5=-P_p*(Grad_g2_q5*Q_p);
Psi4_q5=Psi4_q5(1:2);
Psi4_q6=-P_p*(Grad_g2_q6*Q_p);
Psi4_q6=Psi4_q6(1:2);
Psi4_q7=-P_p*(Grad_g2_q7*Q_p);
Psi4_q7=Psi4_q7(1:2);
Psi4_q8=-P_p*(Grad_g2_q8*Q_p);
Psi4_q8=Psi4_q8(1:2);
Psi4_q9=-P_p*(Grad_g2_q9*Q_p);
Psi4_q9=Psi4_q9(1:2);
Psi4_q10=-P_p*(Grad_g2_q10*Q_p);
Psi4_q10=Psi4_q10(1:2);
Psi4_q11=-P_p*(Grad_g2_q11*Q_p);
Psi4_q11=Psi4_q11(1:2);
Psi4_q12=-P_p*(Grad_g2_q12*Q_p);
Psi4_q12=Psi4_q12(1:2);
Psi4_qp1=Rp_qp1(1:2,1:2)*Param.CM2;
Psi4_qp2=[1;0];
Psi4_qp3=[0;1];


Psi_q=[Psi1_q1 Psi1_q2 Psi1_q3 Psi1_q4 Psi1_q5 Psi1_q6 Psi1_q7 Psi1_q8 Psi1_q9 Psi1_q10 Psi1_q11 Psi1_q12 Psi1_qp1 Psi1_qp2 Psi1_qp3
    Psi2_q1 Psi2_q2 Psi2_q3 Psi2_q4 Psi2_q5 Psi2_q6 Psi2_q7 Psi2_q8 Psi2_q9 Psi2_q10 Psi2_q11 Psi2_q12 Psi2_qp1 Psi2_qp2 Psi2_qp3
    Psi3_q1 Psi3_q2 Psi3_q3 Psi3_q4 Psi3_q5 Psi3_q6 Psi3_q7 Psi3_q8 Psi3_q9 Psi3_q10 Psi3_q11 Psi3_q12 Psi3_qp1 Psi3_qp2 Psi3_qp3
    Psi4_q1 Psi4_q2 Psi4_q3 Psi4_q4 Psi4_q5 Psi4_q6 Psi4_q7 Psi4_q8 Psi4_q9 Psi4_q10 Psi4_q11 Psi4_q12 Psi4_qp1 Psi4_qp2 Psi4_qp3];
res1 =Keps*q - Hq*Param.Forces_Tendons + Psi_q'*lambda + Fe;
res2=Psi;
res=[res1;res2];
end