function ddispdx=ddispdx(m,x,rho,eta,detadx,E,SIGMA,T)
a0=[0.91056314451539 0.63612814494991 2.68613478913903 -26.5473624914884 97.7592087835073 -159.591540865600 91.2977740839123];
a1=[-0.30840169182720 0.18605311591713 -2.50300472586548 21.4197936296668 -65.2558853303492 83.3186804808856 -33.7469229297323];
a2=[-0.09061483509767 0.45278428063920 0.59627007280101 -1.72418291311787 -4.13021125311661 13.7766318697211 -8.67284703679646];
b0=[0.72409469413165 1.11913959304690*2 -1.33419498282114*3 -5.25089420371162*4 5.37112827253230*5 34.4252230677698*6 -50.8003365888685*7];
b1=[-0.57554980753450 0.34975477607218*2 1.29752244631769*3 -4.30386791194303*4 38.5344528930499*5 -26.9710769414608*6 -23.6010990650801*7];
b2=[0.09768831158356 -0.12787874908050*2 -3.05195205099107*3 5.16051899359931*4 -7.76088601041257*5 15.6044623461691*6 -4.23812936930675*7];
n=max(size(m));
mave=0;
for j=1:n
    mave=mave+x(j)*m(j);
end
am=zeros(1,7);
bm=zeros(1,7);
I1=0;
I2=0;
for j=0:6
    am(j+1)=a0(j+1)+((mave-1)/mave)*a1(j+1)+((mave-1)/mave)*((mave-2)/mave)*a2(j+1);
    bm(j+1)=b0(j+1)+((mave-1)/mave)*b1(j+1)+((mave-1)/mave)*((mave-2)/mave)*b2(j+1);
    I1=I1+am(j+1)*(eta^(j));
    I2=I2+bm(j+1)*(eta^(j));
end
damdx=zeros(n,7);
dbmdx=zeros(n,7);
dI1dx=zeros(1,n);
dI2dx=zeros(1,n);
dm2esigma3dx=zeros(1,n);
dm2e2sigma3dx=zeros(1,n);
dC1dx=zeros(1,n);
ddispdx=zeros(1,n);
for k=1:n
    for i=1:7
        damdx(k,i)=m(k)/(mave)^2*a1(i)+m(k)/((mave)^2)*(3-4/mave)*a2(i);
        dbmdx(k,i)=m(k)/(mave)^2*b1(i)+m(k)/((mave)^2)*(3-4/mave)*b2(i);
        dI1dx(1,k)=dI1dx(1,k)+am(i)*(i-1)*detadx(k)*eta^(i-2)+damdx(k,i)*eta^(i-1);
        dI2dx(1,k)=dI2dx(1,k)+bm(i)*(i-1)*detadx(k)*eta^(i-2)+dbmdx(k,i)*eta^(i-1);
    end
    for j=1:n
        dm2esigma3dx(1,k)=dm2esigma3dx(1,k)+2*m(k)*x(j)*m(j)*((E(k,j)/T)^1)*(SIGMA(k,j)^3);
        dm2e2sigma3dx(1,k)=dm2e2sigma3dx(1,k)+2*m(k)*x(j)*m(j)*((E(k,j)/T)^2)*(SIGMA(k,j)^3);
    end
    dC1dx(1,k)=C2(mave,eta)*detadx(k)-(C1(mave,eta)^2)*(m(k)*((8*eta-2*eta^2)/((1-eta)^4))-m(k)*((20*eta-27*eta^2+12*eta^3-2*eta^4)/(((1-eta)*(2-eta))^2)));
    ddispdx(1,k)=-2*pi*rho*(dI1dx(1,k)*m2esigma3(x,m,E,SIGMA,T)+I1*dm2esigma3dx(1,k))-pi*rho*((m(k)*C1(mave,eta)*I2+mave*dC1dx(1,k)*I2+mave*C1(mave,eta)*dI2dx(1,k))*m2e2sigma3(x,m,E,SIGMA,T)+mave*C1(mave,eta)*I2*dm2e2sigma3dx(1,k));
end