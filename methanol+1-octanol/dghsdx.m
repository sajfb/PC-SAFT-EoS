function dghsdx=dghsdx(xi,dxidx,di,dj)
dij=(di*dj)/(di+dj);
dghsdx=(dxidx(3))/(1-xi(3))^2+dij*((3*dxidx(2))/(1-xi(3))^2+(6*xi(2)*dxidx(3))/(1-xi(3))^3)+(dij^2)*((4*xi(2)*dxidx(2))/(1-xi(3))^3+(6*((xi(2))^2)*dxidx(3))/(1-xi(3))^4);
end