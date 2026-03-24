function [res,Grad] = Static_gradient(var,Param)
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
fFe1q=zeros(Param.n*Param.na,Param.n*Param.na);
SFe1q=zeros(Param.n*Param.na,Param.n*Param.na);
fFe2q=zeros(Param.n*Param.na,Param.n*Param.na);
SFe2q=zeros(Param.n*Param.na,Param.n*Param.na);
H1=Ftendons(Y(1),1,q1,Param);
H2=Ftendons(Y(1),1,q2,Param);
fFa10=-phival'*Param.B'*H1;
fFa20=-phival'*Param.B'*H2;
fFe10=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
fFe20=-J2'*Param.M*Ad_function(Inverse_T(g2))*Param.g;
fFe1q0=fFe1q;
fFe2q0=fFe2q;
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

    H1=Ftendons(Y(k),k,q1,Param);
    H2=Ftendons(Y(k),k,q2,Param);
    fFa1=-phival'*Param.B'*H1;
    fFa2=-phival'*Param.B'*H2;
    fFe1=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
    fFe2=-J2'*Param.M*Ad_function(Inverse_T(g2))*Param.g;
    %fFa=-phival'*Param.B'*(Ftendonsval);
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
    P1=N1;
    P2=N2;
    for i=1:Param.n*Param.na
        Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J1(:,i));
        Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))=g2*Hat(J2(:,i));
        Ad_g1_q=[P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R zeros(3);(Skew_symmetric(P_p*(g1*Q_p))*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R))+ (Skew_symmetric(P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p))*(P_R*g1*Q_R)) P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R];
        Ad_g2_q=[P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R zeros(3);(Skew_symmetric(P_p*(g2*Q_p))*(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R))+ (Skew_symmetric(P_p*(Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_p))*(P_R*g2*Q_R)) P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R];
        a1=-(Skew_symmetric((P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'*(P_p*(g1*Q_p))+(P_R*g1*Q_R)'*(P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p)))*(P_R*g1*Q_R)'+Skew_symmetric((P_R*g1*Q_R)'*(P_p*(g1*Q_p)))*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)');
        a2=-(Skew_symmetric((P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)'*(P_p*(g2*Q_p))+(P_R*g2*Q_R)'*(P_p*(Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_p)))*(P_R*g2*Q_R)'+Skew_symmetric((P_R*g2*Q_R)'*(P_p*(g2*Q_p)))*(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)');
        Ad_inv_g1_q=[(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)' zeros(3);a1 (P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'];
        Ad_inv_g2_q=[(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)' zeros(3);a2 (P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)'];

        N1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_g1_q*Param.B*phival;
        I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+(P1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+N1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)))*Param.dX/2;
        J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_inv_g1_q*I1+Ad_function(Inverse_T(g1))*I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i));

        N2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))=Ad_g2_q*Param.B*phival;
        I2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))=I2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))+(P2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))+N2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i)))*Param.dX/2;
        J2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))=Ad_inv_g2_q*I2 + Ad_function(Inverse_T(g2))*I2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i));
        fFe1q(:,i)=-( J1'*Param.M*Ad_inv_g1_q*Param.g + J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))'*Param.M*Ad_function(Inverse_T(g1))*Param.g );
        fFe2q(:,i)=-( J2'*Param.M*Ad_inv_g2_q*Param.g + J2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))'*Param.M*Ad_function(Inverse_T(g2))*Param.g );
        % Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J1(:,i));
        % Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))=g2*Hat(J2(:,i));
        % %Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J1(:,i));
        % %Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))=g2*Hat(J2(:,i));
        % 
        % Ad_g1_q(:,6*(i-1)+1:6*(i))=[P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R zeros(3);(Skew_symmetric(P_p*(g1*Q_p))*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R))+ (Skew_symmetric(P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p))*(P_R*g1*Q_R)) P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R];
        % Ad_g2_q(:,6*(Param.n*Param.na+i-1)+1:6*(Param.n*Param.na+i))=[P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R zeros(3);(Skew_symmetric(P_p*(g2*Q_p))*(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R))+ (Skew_symmetric(P_p*(Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_p))*(P_R*g2*Q_R)) P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R];
        % a1=-(Skew_symmetric((P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'*(P_p*(g1*Q_p))+(P_R*g1*Q_R)'*(P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p)))*(P_R*g1*Q_R)'+Skew_symmetric((P_R*g1*Q_R)'*(P_p*(g1*Q_p)))*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)');
        % a2=-(Skew_symmetric((P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)'*(P_p*(g2*Q_p))+(P_R*g2*Q_R)'*(P_p*(Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_p)))*(P_R*g2*Q_R)'+Skew_symmetric((P_R*g2*Q_R)'*(P_p*(g2*Q_p)))*(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)');
        % Ad_inv_g1_q(:,6*(i-1)+1:6*(i))=[(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)' zeros(3);a1 (P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'];
        % Ad_inv_g2_q(:,6*(Param.n*Param.na+i-1)+1:6*(Param.n*Param.na+i))=[(P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)' zeros(3);a2 (P_R*Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))*Q_R)'];
        % 
        % %P1q=Ad_inv_g1_q(:,6*(i-1)+1:6*(i))*Param.B*phival;
        % N1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_g1_q(:,6*(i-1)+1:6*(i))*Param.B*phival;
        % N2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))=Ad_g2_q(:,6*(Param.n*Param.na+i-1)+1:6*(Param.n*Param.na+i))*Param.B*phival;
        % I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+(P1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+N1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)))*Param.dX/2;
        % I2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))=I2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))+(P2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))+N2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i)))*Param.dX/2;
        % J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_inv_g1_q(:,6*(i-1)+1:6*(i))*I1+Ad_function(Inverse_T(g1))*I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i));
        % J2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i))=Ad_inv_g2_q(:,6*(Param.n*Param.na+i-1)+1:6*(Param.n*Param.na+i))*I2 + Ad_function(Inverse_T(g2))*I2q(:,Param.n*Param.na*(Param.n*Param.na+i-1)+1:Param.n*Param.na*(Param.n*Param.na+i));
    end
    if k==Param.n_X+1
        fFe1qn=fFe1q;
        fFe2qn=fFe2q;
       
    end
    if and(k~=Param.n_X+1,k~=1)
        SFe1q=SFe1q+fFe1q;
        SFe2q=SFe2q+fFe2q;  
    end
    P1q=N1q;
    P2q=N2q;
end
r1=g1(1:3,4);
r2=g2(1:3,4);
R1=g1(1:3,1:3);
R2=g2(1:3,1:3);
Fe1q=Param.L/(Param.n_X)*((fFe1q0+fFe1qn)/2 + SFe1q);
Fe2q=Param.L/(Param.n_X)*((fFe2q0+fFe2qn)/2 + SFe2q);
Feq=[Fe1q, zeros(Param.n*Param.na,Param.n*Param.na+3);zeros(Param.n*Param.na), Fe2q zeros(Param.n*Param.na,3);zeros(3,2*Param.n*Param.na+3)];
Fe1=Param.L/(Param.n_X)*((fFe10+fFe1n)/2 + SFe1);
Fe2=Param.L/(Param.n_X)*((fFe20+fFe2n)/2 + SFe2);

gp=[Rp,rp;0 0 0 1];
%Fep=-Param.Mp*Ad_function(Inverse_T(gp))*Param.g;
%Fep=-Param.Mp*Param.g;
%Fep=[0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0]*Fep;
Fep=[0;-Param.mp*Param.g(4:5)];
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
    Psi1q_q_qp(:,2*Param.n*Param.na*(i-1)+i*3-2:2*Param.n*Param.na*(i)+i*3)=[Psi1_qq(:,2*Param.n*Param.na*(i-1)+1:2*Param.n*Param.na*i) Psi1_qp1q(:,i) Psi1_qp2q(:,i) Psi1_qp3q(:,i)];
    Psi2q_q_qp(:,2*Param.n*Param.na*(i-1)+i*3-2:2*Param.n*Param.na*(i)+i*3)=[Psi2_qq(:,2*Param.n*Param.na*(i-1)+1:2*Param.n*Param.na*(i)) Psi2_qp1q(:,i) Psi2_qp2q(:,i) Psi2_qp3q(:,i)];
    Psi3q_q_qp(:,2*Param.n*Param.na*(i-1)+i*3-2:2*Param.n*Param.na*(i)+i*3)=[Psi3_qq(:,2*Param.n*Param.na*(i-1)+1:2*Param.n*Param.na*(i)) Psi3_qp1q(:,i) Psi3_qp2q(:,i) Psi3_qp3q(:,i)];
    Psi4q_q_qp(:,2*Param.n*Param.na*(i-1)+i*3-2:2*Param.n*Param.na*(i)+i*3)=[Psi4_qq(:,2*Param.n*Param.na*(i-1)+1:2*Param.n*Param.na*(i)) Psi4_qp1q(:,i) Psi4_qp2q(:,i) Psi4_qp3q(:,i)];
end

Psi1_qp1qp1=Hat_inv(Rp_qp1_qp1*Rc1*(P_R*g1*Q_R)'- (P_R*g1*Q_R)*Rc1'*Rp_qp1_qp1');
Psi2_qp1qp1=-Rp_qp1_qp1*[Param.CM1;0];
Psi3_qp1qp1=Hat_inv(Rp_qp1_qp1*Rc2*(P_R*g2*Q_R)'- (P_R*g2*Q_R)*Rc2'*Rp_qp1_qp1');
Psi4_qp1qp1=-Rp_qp1_qp1*[Param.CM2;0];


Psi1q_q_qp=[Psi1q_q_qp Psi1_qqp1 Psi1_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi1_qp1qp2 Psi1_qp1qp3 Psi1_qp2q Psi1_qp2qp1 Psi1_qp2qp2 Psi1_qp2qp3 Psi1_qp3q Psi1_qp3qp1 Psi1_qp3qp2 Psi1_qp3qp3];
Psi1q_q_qp=Psi1q_q_qp(3,:);
Psi2q_q_qp=[Psi2q_q_qp Psi2_qqp1 Psi2_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi2_qp1qp2 Psi2_qp1qp3 Psi2_qp2q Psi2_qp2qp1 Psi2_qp2qp2 Psi2_qp2qp3 Psi2_qp3q Psi2_qp3qp1 Psi2_qp3qp2 Psi2_qp3qp3];
Psi2q_q_qp=Psi2q_q_qp(1:2,:);
Psi3q_q_qp=[Psi3q_q_qp Psi3_qqp1 Psi3_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi3_qp1qp2 Psi3_qp1qp3 Psi3_qp2q Psi3_qp2qp1 Psi3_qp2qp2 Psi3_qp2qp3 Psi3_qp3q Psi3_qp3qp1 Psi3_qp3qp2 Psi3_qp3qp3];
Psi3q_q_qp=Psi3q_q_qp(3,:);
Psi4q_q_qp=[Psi4q_q_qp Psi4_qqp1 Psi4_qp1qp1 zeros(3,8+4*Param.n*Param.na)];%Psi4_qp1qp2 Psi4_qp1qp3 Psi4_qp2q Psi4_qp2qp1 Psi4_qp2qp2 Psi4_qp2qp3 Psi4_qp3q Psi4_qp3qp1 Psi4_qp3qp2 Psi4_qp3qp3];
Psi4q_q_qp=Psi4q_q_qp(1:2,:);

Psi_q_q=[Psi1q_q_qp;Psi2q_q_qp;Psi3q_q_qp;Psi4q_q_qp];
Psi_q_q_transpose_lambda=reshape(lambda'*Psi_q_q,[2*Param.n*Param.na+3, 2*Param.n*Param.na+3]);
Grad= [Keps + Psi_q_q_transpose_lambda + Feq  Psi_q'; Psi_q zeros(6,6)];
res1 =Keps*q - Hq*Param.Forces_Tendons + Psi_q'*lambda + Fe;
res2=Psi;
res=[res1;res2];

% ress=res'*res
% Grad=2*res'*Grad;
end