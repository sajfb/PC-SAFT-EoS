%% 1:CO2 2:MEA 3:H2O 4:HCO3- 5:CO3-2 6:MEACOO- 7:MEAH+ 8:H3O+ 9:OH-
%% reaction1: 2H2O <===> H3O+ + OH-
%% reaction2: CO2 + 2H2O <===> H3O+ + HCO3-
%% reaction3: HCO3- + H2O <===> H3O+ + CO3-2
%% reaction4: MEACOO- + H2O <===> MEA + HCO3-
%% reaction5: MEAH+ + H2O <===> H3O+ + MEA
%%
function [x n]=ChemicalEquilibriumSmithIdeal(n0,T)
% n0=[0.031606090373281 0.319253438113949 4.467258601553829];
% T=313;
nComponents=9;nBalances=4;
R=8.314;
lnKeq1=(132.899-13445.9/T-22.4773*log(T));
lnKeq2=(231.456-12092.1/T-36.7816*log(T));
lnKeq3=(216.049-12431.7/T-35.4819*log(T));
lnKeq4=(2.8898-3635.09/T);
lnKeq5=(2.1211-8189.38/T-0.007484*T);
a=[1 0 0 1 1 1 0 0 0
   0 1 0 0 0 1 1 0 0
   0 0 1 1 1 0 0 1 1
   0 0 0 -1 -2 -1 1 1 -1];
b=[n0(1) n0(2) n0(3) 0]';
bm=b;
Coefficient=[0 0 -2 0 0 0 0 1 1;-1 0 -2 1 0 0 0 1 0;0 0 -1 -1 1 0 0 1 0;0 1 -1 1  0 -1 0 0 0;0 1 -1 0 0 0 -1 1 0;eye(nBalances,nComponents)];
MatlnKx=[lnKeq1 lnKeq2 lnKeq3 lnKeq4 lnKeq5 zeros(1,nBalances)]';
Mu0=-R*T*(Coefficient\MatlnKx);
%n=[0.001168602699100 0.305043956454796 2.660135827212792 0.114557169732994 0.000000981512183 0.114559937469339 0.000000000201929 0.000000804913910]';
n=[n0(1) n0(2) n0(3) 0.114557169732994 0.000000981512183 0.114559937469339 0.114559937469339 0.000000000201929 0.000000804913910]';
Mu=Mu0;
dn=zeros(nComponents,1);
ERROR=1e-13;
ern=zeros(nComponents,1);
er2=100;
er1=1;
while er1>ERROR
    A=zeros(nBalances,nBalances);
    B=zeros(nBalances,1);
    for k=1:nBalances
        for i=1:nBalances
            for j=1:nComponents
                A(k,i)=A(k,i)+a(k,j)*a(i,j)*n(j);
            end
        end
        for j=1:nComponents
            B(k)=B(k)+a(k,j)*n(j)*Mu(j)/R/T;
        end
        B(k)=B(k)+(b(k)-bm(k));
    end
    Psi=A\B;
    for j=1:nComponents
        sumaPsi=0;
        for i=1:nBalances
            sumaPsi=sumaPsi+a(i,j)*Psi(i);
        end
        dn(j)=n(j)*(-Mu(j)/R/T+sumaPsi);
    end
    w=1;
    n1=n;
    n=n+w*dn;
    while min(n)<0
        w=w/2;
        n=n1+w*dn;
    end
    for j=1:nComponents
        Mu(j)=Mu0(j)+R*T*log(n(j)/sum(n));
    end
    for i=1:nComponents
        ern(i)=abs(dn(i)/n(i));
    end
    er1=max(ern);
    er3=abs(er2-er1);
    if er3<ERROR
        break
    end
    er2=er1;
    bm=zeros(nBalances,1);
    for k=1:nBalances
        for i=1:nComponents
            bm(k)=bm(k)+a(k,i)*n(i);
        end
    end
end
x=zeros(nComponents,1);
for i=1:nComponents
    x(i,1)=n(i)/sum(n);
end