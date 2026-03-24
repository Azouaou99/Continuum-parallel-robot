function [q1,f]=Newton_r(q0,P_des,err,Param)
i=0;
f=10;
grad=100*eye(Param.n*Param.na+3+Param.n_lambda);
q1=q0;
while max(abs(f))>err & i<15
   [f,grad]= Static_Gradient_Limit_actuation_V3(q0,P_des,Param);
    q1=q0-(grad)^(-1)*f;
    i=i+1;
    q0=q1;
end
end