function Delta=Delta(di,dj,xi,kapaAB,eAB,T)
dij=(di+dj)/2;
Delta=(dij^3)*ghs(di,dj,xi)*kapaAB*(exp(eAB/T)-1);
end