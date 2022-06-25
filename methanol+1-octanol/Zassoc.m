function Zassoc=Zassoc(nosite,Mole,x,d,xi,rho,X,DELTA)
n=size(Mole,2);
sumZassoc=0;
for i=1:n
    Mole1=cell2mat(Mole(i));
    XAi=cell2mat(X(i));
    for j=1:n
        Mole2=cell2mat(Mole(j));
        XBj=cell2mat(X(j));
        for Ai=1:nosite(i)
            for Bj=1:nosite(j)
                if Mole1(Ai)~=Mole2(Bj)
                    sumZassoc=sumZassoc+x(i)*XAi(Ai)*x(j)*XBj(Bj)*DELTA(i,j)*(1+rhodghsdrho(d(i),d(j),xi)/ghs(d(i),d(j),xi));
                end
            end
        end
    end
end
Zassoc=-rho*sumZassoc/2;