R = [12  13.2 -7.77 -7.77 7     2     4
    2.9 1.45  0.3  -1.19 -3.84 -6.77 4.08
     6.5 5.9  -2.75 -2.75 -2.75 0     -2.75
     3.7 0    -1.5  -1.5  -1.5  0     -1.5];

Fmax = [143.5;127.75;39.2;48.65;12.6;145.6;56];
Fmax=diag(Fmax);

syms l1 l2 l3 th1 th2 th3 ph1;

q=[th1;ph1;th2;th3];
  
r=[l1*cos(ph1)*cos(th1)+l2*cos(ph1)*cos(th1+th2)+l3*cos(ph1)*cos(th1+th2+th3);
   l1*cos(ph1)*sin(th1)+l2*cos(ph1)*sin(th1+th2)+l3*cos(ph1)*sin(th1+th2+th3);
   -l1*sin(ph1)*cos(th1)-l2*sin(ph1)*cos(th1+th2)-l3*sin(ph1)*cos(th1+th2+th3);
   th1+th2+th3+ph1];

J=jacobian(r,q);
transpose(J)


l1=50;l2=31;l3=16;       %meters

th1=pi/4;th2=pi/4;th3=pi/18;ph1=0;  %radians

J=eval(J);
JT=transpose(J);
JTinv=inv(JT);
chi=JTinv*R*Fmax;

fx=chi(1,:);
fy=chi(2,:); %fy is -ve because linprog minimises
fz=chi(3,:);
te=chi(4,:);  

A=[fx;te;fz];
b=[0;0;0];
lb=zeros(7,1);
ub=[1;1;1;1;1;1;1];

%fx 

[x,fval,exitflag,output,lambda]=linprog(-fy,[],[],A,b,lb,ub);
x
fval

//The other objective function
function f=objfun1(x)
syms th1 th2 th3;% a1 a2 a3 a4 a5 a6 a7;

R = [12  13.2 -7.77 -7.77 7     2     4
     6.5 5.9  -2.75 -2.75 -2.75 0     -2.75
     3.7 0    -1.5  -1.5  -1.5  0     -1.5];

Fmax = [143.5;127.75;39.2;48.65;12.6;145.6;56];
Fmax=diag(Fmax);

q=[th1;th2;th3];

r=[ 50*cos(th1)+31*cos(th1+th2)+16*cos(th1+th2+th3);
    50*sin(th1)+31*sin(th1+th2)+16*sin(th1+th2+th3);
    th1+th2+th3];

J=jacobian(r,q);
J=eval(J);
JT=transpose(J);
JTinv=inv(JT);
chi=JTinv*R*Fmax;

f1=chi(1,:);
f2=(-1)*chi(2,:);
t=chi(3,:); 

th1=x(1);
th2=x(2);
th3=x(3);
%a1=x(4);
%a2=x(5);
%a3=x(6);
%a4=x(7);
%a5=x(8);
%a6=x(9);
%a7=x(10);

f=(f2)*[1;1;1;1;1;1;1];
f=eval(f);
%fx=vpa(f1)*[a1;a2;a3;a4;a5;a6;a7];
%te=vpa(t)*[a1;a2;a3;a4;a5;a6;a7];
