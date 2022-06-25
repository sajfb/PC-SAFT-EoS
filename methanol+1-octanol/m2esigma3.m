function m2esigma3=m2esigma3(x,m,E,SIGMA,T)
n=size(x,2);
m2esigma3=0;
for i=1:n
    for j=1:n
        m2esigma3=m2esigma3+x(i)*x(j)*m(i)*m(j)*((E(i,j)/T)^1)*(SIGMA(i,j)^3);
    end
end