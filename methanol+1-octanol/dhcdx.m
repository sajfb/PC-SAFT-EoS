function dahcdx=dhcdx(m,x,xi,dxidx,d)
n=max(size(m));
mave=0;
for i=1:n
    mave=mave+x(i)*m(i);
end
dahsdx=zeros(1,n);
dahcdx=zeros(1,n);
for k=1:n
    dahsdx(1,k)=-dxidx(k,4)/xi(4)*areshs(xi)+1/xi(4)*((3*(dxidx(k,1)*xi(2)+xi(1)*dxidx(k,2)))/(1-xi(3))+(3*xi(1)*xi(2)*dxidx(k,3))/((1-xi(3))^2)+3*(xi(2)^2)*dxidx(k,2)/(xi(3)*((1-xi(3))^2))+(xi(2)^3)*dxidx(k,3)*(3*xi(3)-1)/((xi(3)^2)*((1-xi(3))^3))+((3*((xi(2))^2)*dxidx(k,2)*xi(3)-2*((xi(2))^3)*dxidx(k,3))/(((xi(3))^3))-dxidx(k,4))*log(1-xi(3))+(xi(4)-(((xi(2))^3))/(((xi(3))^2)))*dxidx(k,3)/((1-xi(3))));
    sumdahcdx=0;
    for i=1:n
        sumdahcdx=sumdahcdx+x(i)*(m(i)-1)*dghsdx(xi,dxidx(k,:),d(i),d(i))/ghs(d(i),d(i),xi);
    end
    dahcdx(k)=m(k)*areshs(xi)+mave*dahsdx(1,k)-sumdahcdx-(m(k)-1)*log(ghs(d(k),d(k),xi));
end