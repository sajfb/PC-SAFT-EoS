% DISCLAIMER
% ==========
% 
% You are not allowed to use this code and its associated functions for any
% work which requires verfication (pharmaceutical, medical research, etc.)
%%	This code was written by Sadjad Fakouri for extention PC-SAFT EOS.	%
%%  Feel free to contact with me: sfakor@gmail.com
%%
clear
format long
clc
tic
P=101.32*1000;
Texp=[467.85 461.88 455.06 449.29 444.51 437.95 433.78 424.43 414.68 401.33 388.4 383.64 377.25 368.4 360.18 355.64 352.24 349.62 347.13 344.03 341.99 340.06 338.87 338.15 337.75]';
xMethanol=[0 0.0032 0.0079 0.0127 0.0169 0.0242 0.0278 0.0419 0.0592 0.1002 0.1477 0.1635 0.2064 0.2646 0.3662 0.439 0.5053 0.5823 0.6433 0.7358 0.8143 0.8978 0.9522 0.9819 1]';
yMethanol=[0 0.1564 0.3085 0.4072 0.5053 0.6062 0.6674 0.7691 0.8432 0.9109 0.9512 0.9632 0.9756 0.9848 0.9917 0.9937 0.9957 0.9968 0.9973 0.9979 0.9987 0.999 0.9992 0.9995 1]';
npexp=max(size(Texp));
n=2;
xexp=zeros(npexp,n);
yexp=zeros(npexp,n);
for i=1:npexp
    xexp(i,1)=xMethanol(i);
    xexp(i,2)=1-xMethanol(i);
    yexp(i,1)=yMethanol(i);
    yexp(i,2)=1-yMethanol(i);
end
%%%%%%%%%%%%%%%%% for plot %%%%%%%%%%%%%%%%
np=100+1;
x(:,1)=linspace(0,1,np);
x(:,2)=1-x(:,1);
y=x;
T=460;
Y=[0.5 0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yBubbleTcal=zeros(np,n);
BubbleTcal=zeros(np,1);
philBubbleT=zeros(np,n);
phigBubbleT=zeros(np,n);
KBubbleT=zeros(np,n);
Kx=zeros(np,n);
ZlBubbleT=zeros(np,1);
ZgBubbleT=zeros(np,1);
erBubbleT=zeros(np,1);
eryBubbleT=zeros(np,1);
er=1e-10;
for i=1:np
    %%%%%%%%%% Bubble T %%%%%%%%%%
    T0=1.05*T;
    er0=10;
    while i~=-1
        [philBubbleT(i,:),ZlBubbleT(i)]=PCSAFTMethanolOctanolliquid(T,P,x(i,:));
        [phigBubbleT(i,:),ZgBubbleT(i)]=PCSAFTMethanolOctanolgas(T,P,Y);
        for j=1:n
            KBubbleT(i,j)=philBubbleT(i,j)/phigBubbleT(i,j);
            Kx(i,j)=x(i,j)*KBubbleT(i,j);
        end
        sumKx1=sum(Kx(i,:));
        er1=1;
        while er1>er
            for j=1:n
                Y(j)=x(i,j)*KBubbleT(i,j)/sumKx1;
            end
            [phigBubbleT(i,:),ZgBubbleT(i)]=PCSAFTMethanolOctanolgas(T,P,Y);
            for j=1:n
                KBubbleT(i,j)=philBubbleT(i,j)/phigBubbleT(i,j);
                Kx(i,j)=x(i,j)*KBubbleT(i,j);
            end
            sumKx2=sum(Kx(i,:));
            er1=abs(sumKx2-sumKx1);
            sumKx1=sumKx2;
        end
        er2=abs(sumKx1-1);
        if er2<er
            break
        end
        T1=T-er2/((er0-er2)/(T0-T));
        T0=T;
        T=T1;
        er0=er2;
    end
    BubbleTcal(i,1)=T;
    yBubbleTcal(i,:)=Y;
    %%%%%%%%%% end Bubble T %%%%%%%%%%
end
plot(yBubbleTcal(:,1),BubbleTcal,yexp(:,1),Texp,'o',x(:,1),BubbleTcal,xexp(:,1),Texp,'o')
box off
xlabel('x,y (Methanol)')
ylabel('T/K')
title('P=101.32 kPa & Methanol(2B) 1-Octanol(2B)')
toc