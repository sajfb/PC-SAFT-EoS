function [X1 X2] = NewtonX2B2B(rho,x,DELTA)
XA1=0.1;
XA2=0.1;
XB1=0.1;
XB2=0.1;
n=[XA1 XB1 XA2 XB2]';
tol=10^(-15);
error=1;
F = [XA1+rho*XA1*(x(1)*XB1*DELTA(1,1)+x(2)*XB2*DELTA(1,2))-1
     XB1+rho*XB1*(x(1)*XA1*DELTA(1,1)+x(2)*XA2*DELTA(1,2))-1
     XA2+rho*XA2*(x(1)*XB1*DELTA(2,1)+x(2)*XB2*DELTA(2,2))-1
     XB2+rho*XB2*(x(1)*XA1*DELTA(2,1)+x(2)*XA2*DELTA(2,2))-1];
while error > tol
    J = [rho*(DELTA(1,1)*XB1*x(1)+DELTA(1,2)*XB2*x(2))+1 DELTA(1,1)*XA1*rho*x(1) 0 DELTA(1,2)*XA1*rho*x(2)
         DELTA(1,1)*XB1*rho*x(1) rho*(DELTA(1,1)*XA1*x(1)+DELTA(1,2)*XA2*x(2))+1 DELTA(1,2)*XB1*rho*x(2) 0
         0 DELTA(2,1)*XA2*rho*x(1) rho*(DELTA(2,1)*XB1*x(1)+DELTA(2,2)*XB2*x(2))+1 DELTA(2,2)*XA2*rho*x(2)
         DELTA(2,1)*XB2*rho*x(1) 0 DELTA(2,2)*XB2*rho*x(2) rho*(DELTA(2,1)*XA1*x(1)+DELTA(2,2)*XA2*x(2))+1];
    dn = J\(-F);
    n = n+dn;
    XA1=n(1);XB1=n(2);
    XA2=n(3);XB2=n(4);
    F = [XA1+rho*XA1*(x(1)*XB1*DELTA(1,1)+x(2)*XB2*DELTA(1,2))-1
         XB1+rho*XB1*(x(1)*XA1*DELTA(1,1)+x(2)*XA2*DELTA(1,2))-1
         XA2+rho*XA2*(x(1)*XB1*DELTA(2,1)+x(2)*XB2*DELTA(2,2))-1
         XB2+rho*XB2*(x(1)*XA1*DELTA(2,1)+x(2)*XA2*DELTA(2,2))-1];
     error = max(abs(F));
end
X1=[XA1 XB1];
X2=[XA2 XB2];
%syms x1 x2 XA1 XB1 XA2 XB2 DELTA11 DELTA12 DELTA21 DELTA22 rho
%F=[XA1+rho*XA1*(x1*XB1*DELTA11+x2*XB2*DELTA12)-1
 %  XB1+rho*XB1*(x1*XA1*DELTA11+x2*XA2*DELTA12)-1
%   XA2+rho*XA2*(x1*XB1*DELTA21+x2*XB2*DELTA22)-1
%   XB2+rho*XB2*(x1*XA1*DELTA21+x2*XA2*DELTA22)-1];
%V=[XA1 XB1 XA2 XB2];
%jacobian(F,V)