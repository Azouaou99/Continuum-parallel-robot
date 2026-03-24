function res = Static(var,Param)
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
phival=Phi(Param.na,Param.n,0,Param.L);
Y=0:Param.dX:Param.L;
for k=1:Param.n_X+1
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
    phival=Phi(Param.na,Param.n,Y(k+1),Param.L);
    %g1=g1*expm2(Hat(xival1)*(Param.dX),xival1*(Param.dX));
    %g2=g2*expm2(Hat(xival2)*(Param.dX),xival2*(Param.dX));
    %phival=Phi(Param.na,Param.n,Y(k),Param.L);
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
%Fep=-Param.Mp*Ad_function(Inverse_T(gp))*Param.g;
%Fep=-Param.Mp*Param.g;
%Fep=[0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0]*Fep;
F_p=[0;-Param.mp*Param.g(4:5)];
Param.g(4:6)
Fe=[Fe1;Fe2;F_p];
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
%a=null(Psi_q)'
%null(Psi_q)'
%a*Psi_q'
%[0 0 0 1 1 1;0 0 0 1 1 1;0 0 0 1 1 1;0 0 0 1 1 1;0 0 0 1 1 1;0 0 0 1 1 1;1 1 1 0 0 0;1 1 1 0 0 0;1 1 1 0 0 0;1 1 1 0 0 0;1 1 1 0 0 0;1 1 1 0 0 0;0 0 0 0 0 0;1 0 1 1 0 1;1 1 0 1 1 0]*Psi_q
% Hq*Param.Forces_Tendons
res1 =Keps*q - Hq*Param.Forces_Tendons + Psi_q'*lambda + Fe;
res2=Psi;
res=[res1;res2];

%ress=res'*res
%Grad=2*res'*Grad;

end