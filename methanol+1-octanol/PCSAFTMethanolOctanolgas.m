function [phi, Z]=PCSAFTMethanolOctanolgas(T,P,x)
kB=1.3806503*10^(-23);
n=max(size(x));
m=[1.5255 4.3555];
sigma=[3.2300 3.7145];
e=[188.90 262.74];
d=zeros(1,n);
A=0;B=0;C=0;D=0;
for j=1:n
    d(j)=sigma(j)*(1-0.12*exp(-3*e(j)/T));
    D=D+x(j)*m(j)*(d(j)^0);
    A=A+x(j)*m(j)*(d(j)^1);
    B=B+x(j)*m(j)*(d(j)^2);
    C=C+x(j)*m(j)*(d(j)^3);
end
SIGMA=zeros(n,n);
E=zeros(n,n);
K=[0 0.020
    0.020 0];
for i=1:n
    for j=1:n
        SIGMA(i,j)=(sigma(i)+sigma(j))/2;
        E(i,j)=sqrt(e(i)*e(j))*(1-K(i,j));
    end
end
dxidx=zeros(n,4);
mu=zeros(1,n);
ERROR=10e-15;
drhodeta=6/pi/C;
dxideta(1,1)=A/C;
dxideta(1,2)=B/C;
dxideta(1,3)=1;
dxideta(1,4)=D/C;
eta=10e-10;
er=1;
while er>ERROR
    rhocal=6/pi*eta/C; 
    xi(4)=pi/6*rhocal*D; 
    xi(1)=pi/6*rhocal*A; 
    xi(2)=pi/6*rhocal*B; 
    xi(3)=eta; 
    [zhc, dzhcdeta]=Zhc(x,m,d,xi,dxideta);
    [zdisp, dzdispdeta]=Zdisp(x,m,eta,rhocal,drhodeta,E,SIGMA,T);
    Z=1+zhc+zdisp; 
    Pcal=Z*kB*T*rhocal*(10^(30));
    dZdeta=dzhcdeta+dzdispdeta;
    dPcaldeta=kB*T*(10^30)*(rhocal*dZdeta+Z*drhodeta); 
    eta2=eta+(P-Pcal)/dPcaldeta;
    er=abs(eta-eta2); 
    eta=eta2;
end
for j=1:n
    dxidx(j,1)=pi/6*rhocal*m(j)*d(j)^1;
    dxidx(j,2)=pi/6*rhocal*m(j)*d(j)^2;
    dxidx(j,3)=pi/6*rhocal*m(j)*d(j)^3;
    dxidx(j,4)=pi/6*rhocal*m(j)*d(j)^0;
end
ahc=areshc(m,xi,x,d);
adisp=aresdisp(m,x,rhocal,xi(3),E,SIGMA,T);
ares=ahc+adisp;
dahcdx=dhcdx(m,x,xi,dxidx,d);
dadispdx=ddispdx(m,x,rhocal,xi(3),dxidx(:,3),E,SIGMA,T);
daresdx=dahcdx+dadispdx;
sumaresdx=0;
for j=1:n
    sumaresdx=sumaresdx+x(j)*daresdx(j);
end
for j=1:n
    mu(j)=ares+(Z-1)+daresdx(j)-sumaresdx;
end
lnphi=mu-log(Z);
phi=exp(lnphi);
end