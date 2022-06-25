function m2e2sigma3=m2e2sigma3(x,m,E,SIGMA,T)
n=size(x,2);
m2e2sigma3=0;
for i=1:n
    for j=1:n
        m2e2sigma3=m2e2sigma3+x(i)*x(j)*m(i)*m(j)*((E(i,j)/T)^2)*(SIGMA(i,j)^3);
    end
end