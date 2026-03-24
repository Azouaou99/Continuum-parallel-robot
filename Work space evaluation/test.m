%function [a,x]=test(i)
tic
x=0
for i =1:200
x=x+exp(i);
end
toc
tic
z=1:200;
m=sum(z);
toc
tic
syms i
a=double(symsum(exp(i),i,1,200));
toc