function C1=C1(mave,eta)
C1=(1+mave*((8*eta-2*eta^2)/((1-eta)^4))+(1-mave)*((20*eta-27*eta^2+12*eta^3-2*eta^4)/(((1-eta)*(2-eta))^2)))^(-1);
end