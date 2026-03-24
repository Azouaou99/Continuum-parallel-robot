function [res,Grad] = Static_Gradient_plateform_actuation(var,qp23,Param)
global sig1 
global stab
global sig3
%global sing
q1=var(1:Param.n*Param.na);
q2=var(Param.n*Param.na+1:2*Param.n*Param.na);
qp=[var(2*Param.n*Param.na+1);qp23];
Fep1=var(2*Param.n*Param.na+2:2*Param.n*Param.na+3);
%Rp=[cos(qp(1)) -sin(qp(1)) 0; sin(qp(1)) cos(qp(1)) 0; 0 0 1];
%lambda=var(2*Param.n*Param.na+4 : 2*Param.n*Param.na+3+Param.n_lambda);
Rp=[cos(qp(1)) -sin(qp(1)) 0; sin(qp(1)) cos(qp(1)) 0; 0 0 1];
rp=[qp(2);qp(3);0];
lambda=var(2*Param.n*Param.na+4 : 2*Param.n*Param.na+3+Param.n_lambda);
%S1=var(2*Param.n*Param.na+3+Param.n_lambda+1:2*Param.n*Param.na+3+Param.n_lambda+2);
%S2=var(2*Param.n*Param.na+3+Param.n_lambda+3:2*Param.n*Param.na+3+Param.n_lambda+4);
q=[q1;q2;qp];
SFe1=zeros(Param.n*Param.na,1);
SFe2=zeros(Param.n*Param.na,1);
%SFa1=zeros(Param.n*Param.na,1);
%SFa2=zeros(Param.n*Param.na,1);

P_p=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_p=[0;0;0;1];
P_R=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_R=[1 0 0;0 1 0;0 0 1;0 0 0];

g1=Param.g1;
g2=Param.g2;
J1= zeros(6,Param.n*Param.na);
J2= zeros(6,Param.n*Param.na);
I1= zeros(6,Param.n*Param.na);
I2= zeros(6,Param.n*Param.na);
Y=0:Param.dX:Param.L;
Grad_g1_q=zeros(4,2*4*Param.n*Param.na);
Grad_g2_q=zeros(4,2*4*Param.n*Param.na);
Grad_g1_qq=zeros(4,2*4*Param.n*Param.na*2*Param.n*Param.na);
Grad_g2_qq=zeros(4,2*4*Param.n*Param.na*2*Param.n*Param.na);
Ad_g1_q=zeros(6,2*6*Param.n*Param.na);
Ad_g2_q=zeros(6,2*6*Param.n*Param.na);
Ad_inv_g1_q=zeros(6,2*6*Param.n*Param.na);
Ad_inv_g2_q=zeros(6,2*6*Param.n*Param.na);
phival=Phi(Param.na,Param.n,Y(1),Param.L);
P1=Ad_function(g1)*Param.B*phival;
P2=Ad_function(g2)*Param.B*phival;
P1q=zeros(6,Param.n*Param.na*2*Param.n*Param.na);
P2q=zeros(6,Param.n*Param.na*2*Param.n*Param.na);
N1q=zeros(6,Param.n*Param.na*2*Param.n*Param.na);
N2q=zeros(6,Param.n*Param.na*2*Param.n*Param.na);
I1q=zeros(6,Param.n*Param.na*2*Param.n*Param.na);
I2q=zeros(6,Param.n*Param.na*2*Param.n*Param.na);
J1q=zeros(6,Param.n*Param.na*2*Param.n*Param.na);
J2q=zeros(6,Param.n*Param.na*2*Param.n*Param.na);
H1t=Ftendons(Y(1),1,q1,Param);
H1=H1t(:,2);
H2t=Ftendons(Y(1),1,q2,Param);
H2=H2t(:,1);
fFa10=-phival'*Param.B'*H1;
fFa20=-phival'*Param.B'*H2;
fFe10=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
fFe20=-J2'*Param.M*Ad_function(Inverse_T(g2))*Param.g;
% for k=2:Param.n_X+1
%     xiaval1=Param.xi_a0+phival*q1;
%     xival1=Param.B*xiaval1+Param.B_bar*Param.xi_c;
%     xiaval2=Param.xi_a0+phival*q2;
%     xival2=Param.B*xiaval2+Param.B_bar*Param.xi_c;
%     g1=g1*expm(Hat(xival1)*(Param.dX));
%     g2=g2*expm(Hat(xival2)*(Param.dX));
h=Param.dX;
for k=2:Param.n_X+1
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
    % g1=g1*expm2(Hat(xival1)*(Param.dX),xival1*(Param.dX));
    % g2=g2*expm2(Hat(xival2)*(Param.dX),xival2*(Param.dX));
    N1=Ad_function(g1)*Param.B*phival;
    N2=Ad_function(g2)*Param.B*phival;
    I1=I1+(P1+N1)*Param.dX/2;
    I2=I2+(P2+N2)*Param.dX/2;
    J1=Ad_function(Inverse_T(g1))*I1;
    J2=Ad_function(Inverse_T(g2))*I2;

%     H1t=Ftendons(Y(k),k,q1,Param);
%     H1=H1t(:,2);
%     H2t=Ftendons(Y(k),k,q2,Param);
%     H2=H2t(:,1);
%     fFa1=-phival'*Param.B'*H1;
%     fFa2=-phival'*Param.B'*H2;
    fFe1=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
    fFe2=-J2'*Param.M*Ad_function(Inverse_T(g2))*Param.g;
    %fFa=-phival'*Param.B'*(Ftendonsval);
    if k==Param.n_X+1
        fFe1n=fFe1;
        fFe2n=fFe2;
        %fFa1n=fFa1;
        %fFa2n=fFa2;
    end
    if and(k~=Param.n_X+1,k~=1)
        SFe1=SFe1+fFe1;
        SFe2=SFe2+fFe2;
        %SFa1=SFa1+fFa1;
        %SFa2=SFa2+fFa2;
    end

    P1=N1;
    P2=N2;
    for i=1:Param.n*Param.na
        Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J1(:,i));
        Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))=g2*Hat(J2(:,i));
        %Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J1(:,i));
        %Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))=g2*Hat(J2(:,i));

        Ad_g1_q(:,6*(i-1)+1:6*(i))=[P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R zeros(3);(Skew_symmetric(P_p*(g1*Q_p))*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R))+ (Skew_symmetric(P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p))*(P_R*g1*Q_R)) P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R];
        Ad_g2_q(:,6*(Param.n*Param.na+i-1)+1:6*(Param.n*Param.na+i))=[P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R zeros(3);(Skew_symmetric(P_p*(g2*Q_p))*(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R))+ (Skew_symmetric(P_p*(Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_p))*(P_R*g2*Q_R)) P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R];
        a1=-(Skew_symmetric((P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'*(P_p*(g1*Q_p))+(P_R*g1*Q_R)'*(P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p)))*(P_R*g1*Q_R)'+Skew_symmetric((P_R*g1*Q_R)'*(P_p*(g1*Q_p)))*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)');
        a2=-(Skew_symmetric((P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)'*(P_p*(g2*Q_p))+(P_R*g2*Q_R)'*(P_p*(Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_p)))*(P_R*g2*Q_R)'+Skew_symmetric((P_R*g2*Q_R)'*(P_p*(g2*Q_p)))*(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)');
        Ad_inv_g1_q(:,6*(i-1)+1:6*(i))=[(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)' zeros(3);a1 (P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'];
        Ad_inv_g2_q(:,6*(Param.n*Param.na+i-1)+1:6*(Param.n*Param.na+i))=[(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)' zeros(3);a2 (P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)'];

        %P1q=Ad_inv_g1_q(:,6*(i-1)+1:6*(i))*Param.B*phival;
        N1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_g1_q(:,6*(i-1)+1:6*(i))*Param.B*phival;
        N2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))=Ad_g2_q(:,6*(Param.n*Param.na+i-1)+1:6*(Param.n*Param.na+i))*Param.B*phival;
        I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+(P1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+N1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)))*Param.dX/2;
        I2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))=I2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))+(P2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))+N2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i)))*Param.dX/2;
        J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_inv_g1_q(:,6*(i-1)+1:6*(i))*I1+Ad_function(Inverse_T(g1))*I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i));
        J2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))=Ad_inv_g2_q(:,6*(Param.n*Param.na+i-1)+1:6*(Param.n*Param.na+i))*I2 + Ad_function(Inverse_T(g2))*I2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i));
    end
    P1q=N1q;
    P2q=N2q;
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
Fep=[Fep1;0];
Fe=[Fe1;Fe2;Fep];
% if T(1)>= 0 & T(2) >= 0
%H1q= Param.dX.*((fFa10+fFa1n)./2 + SFa1);
%H2q= Param.dX.*((fFa20+fFa2n)./2 + SFa2);
%Hq=[H1q, zeros(Param.n*Param.na,1);zeros(Param.n*Param.na,1), H2q; zeros(3,2)];
% elseif T(1)< 0 & T(2) >= 0
% H1q= Param.dX.*((fFa10+fFa1n)./2 + SFa1);
% H2q= Param.dX.*((fFa20+fFa2n)./2 + SFa2);
% Hq=[-H1q, zeros(Param.n*Param.na,1);zeros(Param.n*Param.na,1), H2q; zeros(3,2)];
% elseif  T(1)>= 0 & T(2) < 0
% H1q= Param.dX.*((fFa10+fFa1n)./2 + SFa1);
% H2q= Param.dX.*((fFa20+fFa2n)./2 + SFa2);
% Hq=[H1q, zeros(Param.n*Param.na,1);zeros(Param.n*Param.na,1), -H2q; zeros(3,2)];
% elseif T(1)< 0 & T(2) < 0
% H1q= Param.dX.*((fFa10+fFa1n)./2 + SFa1);
% H2q= Param.dX.*((fFa20+fFa2n)./2 + SFa2);
% Hq=[-H1q, zeros(Param.n*Param.na,1);zeros(Param.n*Param.na,1), -H2q; zeros(3,2)];  
% else T(1)== 0 & T(2) == 0
%  Hq=zeros(2*Param.n*Param.na+3,2);
% end
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

Psi1_q=zeros(3,2*Param.n*Param.na);
Psi2_q=zeros(3,2*Param.n*Param.na);
Psi3_q=zeros(3,2*Param.n*Param.na);
Psi4_q=zeros(3,2*Param.n*Param.na);
% for i=1:Param.n*Param.na
%     Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J1(:,i));
%     Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))=g2*Hat(J2(:,i));
%     Psi1_q(:,i)=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'-(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)*Rc1'*Rp');
%     Psi2_q(:,i)=-P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p);
%     Psi3_q(:,Param.n*Param.na+i)=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)'- (P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)*Rc2'*Rp');
%     Psi4_q(:,Param.n*Param.na+i)=-P_p*(Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_p);
% 
% end
% Rp_qp1=[-sin(qp(1)) -cos(qp(1)) 0; cos(qp(1)) -sin(qp(1)) 0; 0 0 0];
% Rp_qp1_qp1=[-cos(qp(1)) sin(qp(1)) 0; -sin(qp(1)) -cos(qp(1)) 0; 0 0 0];
% Psi1_q=Psi1_q(3,:);
% Psi2_q=Psi2_q(1:2,:);
% Psi3_q=Psi3_q(3,:);
% Psi4_q=Psi4_q(1:2,:);
% 
% Psi1_qp1=Hat_inv(Rp_qp1*Rc1*(P_R*g1*Q_R)'- (P_R*g1*Q_R)*Rc1'*Rp_qp1');
% Psi1_qp1=Psi1_qp1(3);
% Psi1_qp2=0;
% Psi1_qp3=0;
% 
% Psi2_qp1=-Rp_qp1(1:2,1:2)*Param.CM1;
% Psi2_qp2=[1;0];
% Psi2_qp3=[0;1];
% 
% Psi3_qp1=Hat_inv(Rp_qp1*Rc2*(P_R*g2*Q_R)'- (P_R*g2*Q_R)*Rc2'*Rp_qp1');
% Psi3_qp1=Psi3_qp1(3);
% Psi3_qp2=0;
% Psi3_qp3=0;
% 
% Psi4_qp1=-Rp_qp1(1:2,1:2)*Param.CM2;
% Psi4_qp2=[1;0];
% Psi4_qp3=[0;1];


% Psi_q=[Psi1_q1 Psi1_q2 Psi1_q3 Psi1_q4 Psi1_q5 Psi1_q6 Psi1_q7 Psi1_q8 Psi1_q9 Psi1_q10 Psi1_q11 Psi1_q12 Psi1_qp1 Psi1_qp2 Psi1_qp3
%     Psi2_q1 Psi2_q2 Psi2_q3 Psi2_q4 Psi2_q5 Psi2_q6 Psi2_q7 Psi2_q8 Psi2_q9 Psi2_q10 Psi2_q11 Psi2_q12 Psi2_qp1 Psi2_qp2 Psi2_qp3
%     Psi3_q1 Psi3_q2 Psi3_q3 Psi3_q4 Psi3_q5 Psi3_q6 Psi3_q7 Psi3_q8 Psi3_q9 Psi3_q10 Psi3_q11 Psi3_q12 Psi3_qp1 Psi3_qp2 Psi3_qp3
%     Psi4_q1 Psi4_q2 Psi4_q3 Psi4_q4 Psi4_q5 Psi4_q6 Psi4_q7 Psi4_q8 Psi4_q9 Psi4_q10 Psi4_q11 Psi4_q12 Psi4_qp1 Psi4_qp2 Psi4_qp3];

% Psi_q=[Psi1_q Psi1_qp1 Psi1_qp2 Psi1_qp3
%     Psi2_q Psi2_qp1 Psi2_qp2 Psi2_qp3
%     Psi3_q Psi3_qp1 Psi3_qp2 Psi3_qp3
%     Psi4_q Psi4_qp1 Psi4_qp2 Psi4_qp3];

Psi1_qq=zeros(3,2*Param.n*Param.na*2*Param.n*Param.na);
Psi2_qq=zeros(3,2*Param.n*Param.na*2*Param.n*Param.na);
Psi3_qq=zeros(3,2*Param.n*Param.na*2*Param.n*Param.na);
Psi4_qq=zeros(3,2*Param.n*Param.na*2*Param.n*Param.na);

Psi1_qqp1=zeros(3,2*Param.n*Param.na);
Psi1_qqp2=zeros(3,2*Param.n*Param.na);
Psi1_qqp3=zeros(3,2*Param.n*Param.na);

Psi2_qqp1=zeros(3,2*Param.n*Param.na);
Psi2_qqp2=zeros(3,2*Param.n*Param.na);
Psi2_qqp3=zeros(3,2*Param.n*Param.na);

Psi3_qqp1=zeros(3,2*Param.n*Param.na);
Psi3_qqp2=zeros(3,2*Param.n*Param.na);
Psi3_qqp3=zeros(3,2*Param.n*Param.na);

Psi4_qqp1=zeros(3,2*Param.n*Param.na);
Psi4_qqp2=zeros(3,2*Param.n*Param.na);
Psi4_qqp3=zeros(3,2*Param.n*Param.na);

Psi1_qp1q=zeros(3,2*Param.n*Param.na);
Psi1_qp2q=zeros(3,2*Param.n*Param.na);
Psi1_qp3q=zeros(3,2*Param.n*Param.na);

Psi2_qp1q=zeros(3,2*Param.n*Param.na);
Psi2_qp2q=zeros(3,2*Param.n*Param.na);
Psi2_qp3q=zeros(3,2*Param.n*Param.na);

Psi3_qp1q=zeros(3,2*Param.n*Param.na);
Psi3_qp2q=zeros(3,2*Param.n*Param.na);
Psi3_qp3q=zeros(3,2*Param.n*Param.na);

Psi4_qp1q=zeros(3,2*Param.n*Param.na);
Psi4_qp2q=zeros(3,2*Param.n*Param.na);
Psi4_qp3q=zeros(3,2*Param.n*Param.na);

Rp_qp1=[-sin(qp(1)) -cos(qp(1)) 0; cos(qp(1)) -sin(qp(1)) 0; 0 0 0];
Rp_qp1_qp1=[-cos(qp(1)) sin(qp(1)) 0; -sin(qp(1)) -cos(qp(1)) 0; 0 0 0];
for i=1:Param.n*Param.na
    %Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J1(:,i));
    %Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))=g2*Hat(J2(:,i));
    Psi1_q(:,i)=Hat_inv(Rp*Rc1*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'-(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)*Rc1'*Rp');
    Psi2_q(:,i)=-P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p);
    Psi3_q(:,Param.n*Param.na+i)=Hat_inv(Rp*Rc2*(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)'- (P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)*Rc2'*Rp');
    Psi4_q(:,Param.n*Param.na+i)=-P_p*(Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_p);


    Psi1_qqp1(:,i)=Hat_inv(Rp_qp1*Rc1*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'-(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)*Rc1'*Rp_qp1');
    Psi3_qqp1(:,Param.n*Param.na+i)=Hat_inv(Rp_qp1*Rc2*(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)'- (P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)*Rc2'*Rp_qp1');
    Psi1_qp1q=Psi1_qqp1;
    Psi3_qp1q=Psi3_qqp1;

    for j=1:Param.n*Param.na
        Grad_g1_qq(:,(i-1)*4*Param.n*Param.na+4*(j-1)+1:(i-1)*4*Param.n*Param.na+4*(j))=Grad_g1_q(:,4*(i-1)+1:4*(i))*Hat(J1(:,j))+ g1*Hat(J1q(:,Param.n*Param.na*(i-1)+j));
        Grad_g2_qq(:,(i-1+Param.n*Param.na)*4*Param.n*Param.na+4*(Param.n*Param.na+j-1)+1:(i-1+Param.n*Param.na)*4*Param.n*Param.na+4*(Param.n*Param.na+j))=Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Hat(J2(:,j))+ g2*Hat(J2q(:,Param.n*Param.na*(Param.n*Param.na + i-1)+j));
        Psi1_qq(:,2*Param.n*Param.na*(i-1)+j)=Hat_inv(Rp*Rc1*(P_R*Grad_g1_qq(:,(i-1)*4*Param.n*Param.na+4*(j-1)+1:(i-1)*4*Param.n*Param.na+4*(j))*Q_R)'-(P_R*Grad_g1_qq(:,(i-1)*4*Param.n*Param.na+4*(j-1)+1:(i-1)*4*Param.n*Param.na+4*(j))*Q_R)*Rc1'*Rp');
        Psi2_qq(:,2*Param.n*Param.na*(i-1)+j)=-P_p*(Grad_g1_qq(:,(i-1)*4*Param.n*Param.na+4*(j-1)+1:(i-1)*4*Param.n*Param.na+4*(j))*Q_p);
        Psi3_qq(:,2*Param.n*Param.na*(i-1+Param.n*Param.na)+Param.n*Param.na+j)=Hat_inv(Rp*Rc2*(P_R*Grad_g2_qq(:,(i-1+Param.n*Param.na)*4*Param.n*Param.na+4*(Param.n*Param.na+j-1)+1:(i-1+Param.n*Param.na)*4*Param.n*Param.na+4*(Param.n*Param.na+j))*Q_R)'- (P_R*Grad_g2_qq(:,(i-1+Param.n*Param.na)*4*Param.n*Param.na+4*(Param.n*Param.na+j-1)+1:(i-1+Param.n*Param.na)*4*Param.n*Param.na+4*(Param.n*Param.na+j))*Q_R)*Rc2'*Rp');
        Psi4_qq(:,2*Param.n*Param.na*(i-1+Param.n*Param.na)+Param.n*Param.na+j)=-P_p*(Grad_g2_qq(:,(i-1+Param.n*Param.na)*4*Param.n*Param.na+4*(Param.n*Param.na+j-1)+1:(i-1+Param.n*Param.na)*4*Param.n*Param.na+4*(Param.n*Param.na+j))*Q_p);
    end
end
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

Psi_q=[Psi1_q Psi1_qp1 Psi1_qp2 Psi1_qp3
    Psi2_q Psi2_qp1 Psi2_qp2 Psi2_qp3
    Psi3_q Psi3_qp1 Psi3_qp2 Psi3_qp3
    Psi4_q Psi4_qp1 Psi4_qp2 Psi4_qp3];

for i=1:2*Param.n*Param.na
    %Psi1_qq(:,2*Param.n*Param.na*(i-1)+1:2*Param.n*Param.na*i)
    Psi1q_q_qp(:,2*Param.n*Param.na*(i-1)+i*3-2:2*Param.n*Param.na*(i)+i*3)=[Psi1_qq(:,2*Param.n*Param.na*(i-1)+1:2*Param.n*Param.na*i) Psi1_qp1q(:,i) Psi1_qp2q(:,i) Psi1_qp3q(:,i)];
    Psi2q_q_qp(:,2*Param.n*Param.na*(i-1)+i*3-2:2*Param.n*Param.na*(i)+i*3)=[Psi2_qq(:,2*Param.n*Param.na*(i-1)+1:2*Param.n*Param.na*(i)) Psi2_qp1q(:,i) Psi2_qp2q(:,i) Psi2_qp3q(:,i)];
    Psi3q_q_qp(:,2*Param.n*Param.na*(i-1)+i*3-2:2*Param.n*Param.na*(i)+i*3)=[Psi3_qq(:,2*Param.n*Param.na*(i-1)+1:2*Param.n*Param.na*(i)) Psi3_qp1q(:,i) Psi3_qp2q(:,i) Psi3_qp3q(:,i)];
    Psi4q_q_qp(:,2*Param.n*Param.na*(i-1)+i*3-2:2*Param.n*Param.na*(i)+i*3)=[Psi4_qq(:,2*Param.n*Param.na*(i-1)+1:2*Param.n*Param.na*(i)) Psi4_qp1q(:,i) Psi4_qp2q(:,i) Psi4_qp3q(:,i)];
end

Psi1_qp1qp1=Hat_inv(Rp_qp1_qp1*Rc1*(P_R*g1*Q_R)'- (P_R*g1*Q_R)*Rc1'*Rp_qp1_qp1');
Psi2_qp1qp1=-Rp_qp1_qp1*[Param.CM1;0];
Psi3_qp1qp1=Hat_inv(Rp_qp1_qp1*Rc2*(P_R*g2*Q_R)'- (P_R*g2*Q_R)*Rc2'*Rp_qp1_qp1');
Psi4_qp1qp1=-Rp_qp1_qp1*[Param.CM2;0];


% Psi1q_q_qp=[Psi1q_q_qp Psi1_qqp1 Psi1_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi1_qp1qp2 Psi1_qp1qp3 Psi1_qp2q Psi1_qp2qp1 Psi1_qp2qp2 Psi1_qp2qp3 Psi1_qp3q Psi1_qp3qp1 Psi1_qp3qp2 Psi1_qp3qp3];
% Psi1q_q_qp=Psi1q_q_qp(3,:);
% Psi2q_q_qp=[Psi2q_q_qp Psi2_qqp1 Psi2_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi2_qp1qp2 Psi2_qp1qp3 Psi2_qp2q Psi2_qp2qp1 Psi2_qp2qp2 Psi2_qp2qp3 Psi2_qp3q Psi2_qp3qp1 Psi2_qp3qp2 Psi2_qp3qp3];
% Psi2q_q_qp=Psi2q_q_qp(1:2,:);
% Psi3q_q_qp=[Psi3q_q_qp Psi3_qqp1 Psi3_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi3_qp1qp2 Psi3_qp1qp3 Psi3_qp2q Psi3_qp2qp1 Psi3_qp2qp2 Psi3_qp2qp3 Psi3_qp3q Psi3_qp3qp1 Psi3_qp3qp2 Psi3_qp3qp3];
% Psi3q_q_qp=Psi3q_q_qp(3,:);
% Psi4q_q_qp=[Psi4q_q_qp Psi4_qqp1 Psi4_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi4_qp1qp2 Psi4_qp1qp3 Psi4_qp2q Psi4_qp2qp1 Psi4_qp2qp2 Psi4_qp2qp3 Psi4_qp3q Psi4_qp3qp1 Psi4_qp3qp2 Psi4_qp3qp3];
% Psi4q_q_qp=Psi4q_q_qp(1:2,:);

GradPsi1_q=[Psi1q_q_qp Psi1_qqp1 Psi1_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi1_qp1qp2 Psi1_qp1qp3 Psi1_qp2q Psi1_qp2qp1 Psi1_qp2qp2 Psi1_qp2qp3 Psi1_qp3q Psi1_qp3qp1 Psi1_qp3qp2 Psi1_qp3qp3];
GradPsi1_q=GradPsi1_q(3,:);
GradPsi2_q=[Psi2q_q_qp Psi2_qqp1 Psi2_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi2_qp1qp2 Psi2_qp1qp3 Psi2_qp2q Psi2_qp2qp1 Psi2_qp2qp2 Psi2_qp2qp3 Psi2_qp3q Psi2_qp3qp1 Psi2_qp3qp2 Psi2_qp3qp3];
GradPsi2_q=GradPsi2_q(1:2,:);
GradPsi3_q=[Psi3q_q_qp Psi3_qqp1 Psi3_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi3_qp1qp2 Psi3_qp1qp3 Psi3_qp2q Psi3_qp2qp1 Psi3_qp2qp2 Psi3_qp2qp3 Psi3_qp3q Psi3_qp3qp1 Psi3_qp3qp2 Psi3_qp3qp3];
GradPsi3_q=GradPsi3_q(3,:);
GradPsi4_q=[Psi4q_q_qp Psi4_qqp1 Psi4_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi4_qp1qp2 Psi4_qp1qp3 Psi4_qp2q Psi4_qp2qp1 Psi4_qp2qp2 Psi4_qp2qp3 Psi4_qp3q Psi4_qp3qp1 Psi4_qp3qp2 Psi4_qp3qp3];
GradPsi4_q=GradPsi4_q(1:2,:);

%Psi_q_q=[Psi1q_q_qp;Psi2q_q_qp;Psi3q_q_qp;Psi4q_q_qp];
Psi_q_q=[GradPsi1_q;GradPsi2_q;GradPsi3_q;GradPsi4_q];
Psi_q_q_transpose_lambda=reshape(lambda'*Psi_q_q,[2*Param.n*Param.na+3, 2*Param.n*Param.na+3]);
grad_Fep1=[zeros(2*Param.n*Param.na,2);eye(2);0 0];
grad_Fep=[zeros(2*Param.n*Param.na+3,2*Param.n*Param.na+1) , grad_Fep1];
Psi_q_T=[Psi1_q Psi1_qp1 0 0
    Psi2_q Psi2_qp1 zeros(2)
    Psi3_q Psi3_qp1 0 0
    Psi4_q Psi4_qp1 zeros(2)];
%a=[zeros(4,2*Param.n*Param.na+1) [eye(2);eye(2)] zeros(4,6) [2*S1 zeros(2,1);zeros(2,1) 2*S2]]
%b=[Keps + Psi_q_q_transpose_lambda + grad_H Psi_q' zeros(2*Param.n*Param.na+3,4)]
%Grad= [Keps + Psi_q_q_transpose_lambda + grad_H Psi_q' zeros(2*Param.n*Param.na+3,4);zeros(4,2*Param.n*Param.na+1) [eye(2);eye(2)] zeros(4,6) diag([-sign(S1(1)) -sign(S1(2)) sign(S2(1)) sign(S2(2))]);Psi_q_T zeros(6,10)];%diag([-sign(S1(1)) -sign(S1(2)) sign(S2(1)) sign(S2(2))])] ;
Grad= [[Keps + Psi_q_q_transpose_lambda + grad_Fep , Psi_q']; [Psi_q_T,zeros(6,6)]];%diag([-sign(S1(1)) -sign(S1(2)) sign(S2(1)) sign(S2(2))])] ;
%Grad= [Keps + grad_H Psi_q';Psi_q_T zeros(6,6)];%diag([-sign(S1(1)) -sign(S1(2)) sign(S2(1)) sign(S2(2))])] ;
%Grad= [Keps + Psi_q_q_transpose_lambda + grad_H Psi_q' zeros(2*Param.n*Param.na+3,4); Psi_q_T zeros(6,10);zeros(4,2*Param.n*Param.na+1) [eye(2);eye(2)] zeros(4,6) diag([-2*S1(1) -2*S1(2) 2*S2(1) 2*S2(2)])];%diag([-sign(S1(1)) -sign(S1(2)) sign(S2(1)) sign(S2(2))])] ;
res1 =Keps*q + Psi_q'*lambda + Fe;
%res2=T-Param.T_min - abs(S1);
%res3=T-Param.T_max + abs(S2);
res2=Psi;

%res3=T-Param.T_min - S1.^2;
%res4=T-Param.T_max + S2.^2;
res=[res1;res2];%;res3;res4];

% T=((Param.T_min + abs(S1)) + (Param.T_max - abs(S2)))./2;
% %T=((Param.T_min + S1.^2) + (Param.T_max - S2.^2))./2;
% %T=((Param.T_min + S1) + (Param.T_max - S2))./2;
% res1 =Keps*q - Hq*T + Psi_q'*lambda + Fe;
% res2=Psi;
% res3=(Param.T_min + abs(S1)) - (Param.T_max - abs(S2));
% res=[res1;res2;res3];
eps_ql1=Keps(:,1:2*Param.n*Param.na)+reshape(lambda'*[Psi1q_q_qp(3,:);Psi2q_q_qp(1:2,:);Psi3q_q_qp(3,:);Psi4q_q_qp(1:2,:)],[2*Param.n*Param.na+3, 2*Param.n*Param.na]);
eps_ql=eps_ql1;
%eps_ql=pinv(Keps(:,1:2*Param.n*Param.na))*eps_ql
for i=1:min(size(eps_ql1))
    eps_ql(i,i)=eps_ql1(i,i)/Keps(i,i);
% % if S11(i,i)<10^(-5)
% %     sing=1;
 end
eps_qp=Keps(:,2*Param.n*Param.na+1:2*Param.n*Param.na+3)+reshape(lambda'*[GradPsi1_q(:,length(Psi1q_q_qp(1,:))+1:end);GradPsi2_q(:,length(Psi1q_q_qp(1,:))+1:end);GradPsi3_q(:,length(Psi1q_q_qp(1,:))+1:end);GradPsi4_q(:,length(Psi1q_q_qp(1,:))+1:end)],[2*Param.n*Param.na+3, 3]);
eps_Fep=grad_Fep1;
eps_lambda=Psi_q';
eps_S=zeros(2*Param.n*Param.na+3,4);
Psi_ql=[Psi1_q 
    Psi2_q 
    Psi3_q 
    Psi4_q];
Psi_qp=[Psi1_qp1 Psi1_qp2 Psi1_qp3
     Psi2_qp1 Psi2_qp2 Psi2_qp3
    Psi3_qp1 Psi3_qp2 Psi3_qp3
     Psi4_qp1 Psi4_qp2 Psi4_qp3];
Psi_Fep=zeros(6,2);
%eps_T=-Hq;
Psi_lambda=zeros(6,6);
Psi_S=zeros(6,4);

%S1=((T-Param.T_min)-[0.5;0.5]).^(20);
%S2=((Param.T_max-T-[0.5;0.5])).^(20);
%S1=((T-Param.T_min)/3).^(10);
%S2=((Param.T_max-T)/3).^(10);
% l=2;
% S1=-1/l*log(Param.T_max./(T-Param.T_min)-1)
% S2=-1/l*log(Param.T_max./(Param.T_max-T)-1)
% S_ql=zeros(4,2*Param.n*Param.na);
% S_qp=zeros(4,3);
% S_T=1/3*[eye(2);eye(2)];
% % S_S=diag([-S1(1) -S1(2) S2(1) S2(2)]);
% S_S11=-Param.T_max(1)*(l*exp(-l*S1(1))/(1+exp(-l*S1(1)))^2)
% S_S12=-Param.T_max(2)*(l*exp(-l*S1(2))/(1+exp(-l*S1(2)))^2)
% S_S21=Param.T_max(1)*(l*exp(-l*S2(1))/(1+exp(-l*S2(1)))^2)
% S_S22=Param.T_max(2)*(l*exp(-l*S2(2))/(1+exp(-l*S2(2)))^2)
% S_S=diag([S_S11 S_S12 S_S21 S_S22])
%n_S_S=null(S_S');
%Lambda_S=[Psi_q' zeros(2*Param.n*Param.na+3,4);zeros(4,6) S_S];
%Z=null(Lambda_S');
Z=null(Psi_q);
%Z=[Z(:,1:6)];
%A=[Z Psi_q'];
%r1=det(A);
%Psi_ql
%Z'*eps_ql
% L=[Z'*[eps_ql;S_ql];Psi_ql];
% P=[Z'*[eps_qp;S_qp];Psi_qp];
% M=[Z'*[eps_T;S_T];Psi_T];

%AA=Z'*eps_ql;
%rank(AA,10^(-10))
%n_S_S'*S_T
 L=[Z'*eps_ql1;Psi_ql];
 P=[Z'*eps_qp;Psi_qp];
 M=[Z'*eps_Fep ; Psi_Fep];
 %MS1=[M zeros(length(M(:,1)),4) L;S_T S_S S_ql]
% L=[Z'*eps_ql;Psi_ql];
% P=[Z'*eps_qp;Psi_qp];
% M=[Z'*eps_T;Psi_T];
H=[eps_ql1 eps_qp];
Hr=Z'*H*Z;
stab=min(eig(Hr));


%[U1,S11,V1]=svd(MS1);
%[U2,S12,V2]=svd([L]);
[U1,S11,V1]=svd([P L]);
sig1=cond([P L])^(-1);
%[U4,S14,V4]=svd([M]); 
% sigma=1;

 %if  T(1)<1 | T(2)<1 |T(1)<1 | T(2)>7.5 | T(1)>7.5 %1e-5
%sig1=0;
% else
% for i=1:min(size(S11))
%     sigma1(i)=S11(i,i);
% end
%  sig1=min(sigma1);
% end

% sig2=min(sigma2);
% sig3=min(sigma3);
%invr=T1'*((T1*T1')^(-1))
%r2=rank(S_S,1e-1);
%r2=rank(T1,1e-1)

