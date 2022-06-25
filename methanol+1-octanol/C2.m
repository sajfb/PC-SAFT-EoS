function C2=C2(mave,eta)
C2=-(C1(mave,eta)^2)*(mave*((-4*eta^2+20*eta+8)/((1-eta)^5))+(1-mave)*((2*eta^3+12*eta^2-48*eta+40)/(((1-eta)*(2-eta))^3)));
end