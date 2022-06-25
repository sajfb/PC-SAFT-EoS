function dassocdx=dassocdx(nosite,Mole,x,d,xi,rho,X,DELTA,dxidx)
n=size(x,2);
dassocdx=zeros(1,n);
for k=1:n
    XAk=cell2mat(X(k));
    sumlnXAk=0;
    for Ak=1:nosite(k)
        sumlnXAk=sumlnXAk+log(XAk(Ak));
    end
    sumdassocdx=0;
    for i=1:n
        Mole1=cell2mat(Mole(i));
        XAi=cell2mat(X(i));
        for j=1:n
            Mole2=cell2mat(Mole(j));
            XBj=cell2mat(X(j));
            for Ai=1:nosite(i)
                for Bj=1:nosite(j)
                    if Mole1(Ai)~=Mole2(Bj)
                        sumdassocdx=sumdassocdx+x(i)*XAi(Ai)*x(j)*XBj(Bj)*DELTA(i,j)*(dghsdx(xi,dxidx(k,:),d(i),d(j))/ghs(d(i),d(j),xi));
                    end
                end
            end
        end
    end
    dassocdx(1,k)=sumlnXAk-rho*sumdassocdx/2;
end