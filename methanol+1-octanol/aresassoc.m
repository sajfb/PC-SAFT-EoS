function aassoc=aresassoc(nosite,x,X)
n=size(x,2);
aassoc=0;
for i=1:n
    XAi=cell2mat(X(i));
    for Ai=1:nosite(i)
        aassoc=aassoc+x(i)*(log(XAi(Ai))-XAi(Ai)/2+1/2);
    end
end