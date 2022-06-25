function ahc=areshc(m,xi,x,d)
n=max(size(m));
ahs=areshs(xi);
sumahc=0;
mave=0;
for j=1:n
    sumahc=sumahc+x(j)*(m(j)-1)*log(ghs(d(j),d(j),xi));
    mave=mave+x(j)*m(j);
end
ahc=mave*ahs-sumahc;
end