function [phi, Z]=PCSAFTMethanolOctanolliquid(T,P,x)
kB=1.3806503*1e-23;
n=max(size(x));
m=[1.5255 4.3555];
sigma=[3.2300 3.7145];
e=[188.90 262.74];
kapaAB=[0.035176 0.002197];
eAB=[2899.5 2754.8];
Molecule1=[+1 -1];
Molecule2=[+1 -1];
Mole={Molecule1,Molecule2};
nosite=[size(Molecule1,2) size(Molecule2,2)];
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
KAPA=zeros(n,n);
EAB=zeros(n,n);
K=[0 0.020
    0.020 0];
for i=1:n
    for j=1:n
        if i==j
            SIGMA(i,j)=sigma(i);
            E(i,j)=e(i);
            KAPA(i,j)=kapaAB(i);
            EAB(i,j)=eAB(i);
        else
            SIGMA(i,j)=(sigma(i)+sigma(j))/2;
            E(i,j)=sqrt(e(i)*e(j))*(1-K(i,j));
            KAPA(i,j)=sqrt(kapaAB(i)*kapaAB(j))*(((sqrt(sigma(i)*sigma(j)))/((sigma(i)+sigma(j))/2))^3);
            EAB(i,j)=(eAB(i)+eAB(j))/2;
        end
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
eta=0.4;
er=1;
DELTA=zeros(n,n);
while er>ERROR
    rhocal=6/pi*eta/C;
    xi(4)=pi/6*rhocal*D;
    xi(1)=pi/6*rhocal*A;
    xi(2)=pi/6*rhocal*B;
    xi(3)=eta;
    [zhc, dzhcdeta]=Zhc(x,m,d,xi,dxideta);
    [zdisp, dzdispdeta]=Zdisp(x,m,eta,rhocal,drhodeta,E,SIGMA,T);
    %%%%%%%%%%%%%%%% Z association %%%%%%%%%%%%%%%%%
    for g=1:n
        for h=1:n
            DELTA(g,h)=Delta(d(g),d(h),xi,KAPA(g,h),EAB(g,h),T);
        end
    end
    [X1, X2]=NewtonX2B2B(rhocal,x,DELTA);
    X={X1,X2};
    zassoc=Zassoc(nosite,Mole,x,d,xi,rhocal,X,DELTA);
    %%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%
    Z=1+zhc+zdisp+zassoc;
    Pcal=Z*kB*T*rhocal*(1e30);
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
aassoc=aresassoc(nosite,x,X);
ares=ahc+adisp+aassoc;
dahcdx=dhcdx(m,x,xi,dxidx,d);
dadispdx=ddispdx(m,x,rhocal,xi(3),dxidx(:,3),E,SIGMA,T);
daassocdx=dassocdx(nosite,Mole,x,d,xi,rhocal,X,DELTA,dxidx);
daresdx=dahcdx+dadispdx+daassocdx;
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