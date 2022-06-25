function [Zhc dZhcdeta]=Zhc(x,m,d,xi,dxideta)
n=max(size(x));
Zhs=xi(3)/(1-xi(3))+3*xi(1)*xi(2)/xi(4)/((1-xi(3))^2)+(3*((xi(2))^3)-xi(3)*((xi(2))^3))/xi(4)/((1-xi(3))^3);
mave=0;
sumZhc=0;
sumdZhcdeta=0; 
for j=1:n
    mave=mave+x(j)*m(j);
    sumZhc=sumZhc+x(j)*(m(j)-1)*rhodghsdrho(d(j),d(j),xi)/ghs(d(j),d(j),xi);
    KPrime0ij=(1+xi(3))/((1-xi(3))^3)+3*(d(j)/2)*(2*xi(2)*(xi(3)+2)+dxideta(2)*(1-(xi(3)^2)))/((1-xi(3))^4)+2*(d(j)/2)^2*xi(2)*(3*xi(2)*(xi(3)+3)-2*dxideta(2)*((xi(3)^2)+xi(3)-2))/((1-xi(3))^5);
    K0ij=rhodghsdrho(d(j),d(j),xi);
    dghsdeta=(1/(1-xi(3)))/(1-xi(3))+3*(d(j)/2)*xi(2)*((1+xi(3))/xi(3))/(1-xi(3))^3+(((d(j)/2))*(xi(2)/(1-xi(3)))/(1-xi(3)))^2*(4+2*xi(3))/xi(3);
    dKiideta=KPrime0ij/ghs(d(j),d(j),xi)-K0ij/((ghs(d(j),d(j),xi))^2)*dghsdeta;
    sumdZhcdeta=sumdZhcdeta+x(j)*(m(j)-1)*dKiideta;
end
Zhc=mave*Zhs-sumZhc;
dZhsdeta=1/((1-xi(3))^2)+3*(dxideta(1)*xi(2)*xi(4)+xi(1)*dxideta(2)*xi(4)-xi(1)*xi(2)*dxideta(4))/((xi(4)*(1-xi(3)))^2)+(xi(2)*xi(4)*(6*xi(1)-xi(2)^2)+(3-xi(3))*(3*dxideta(2)*xi(4)-xi(2)*dxideta(4))*(xi(2)^2))/((xi(4)^2)*((1-xi(3))^3))+3*(xi(2)^3)*(3-xi(3))/(xi(4)*((1-xi(3))^4));
dZhcdeta=mave*dZhsdeta-sumdZhcdeta;