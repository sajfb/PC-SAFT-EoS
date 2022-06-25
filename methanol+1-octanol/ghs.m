function ghs=ghs(di,dj,xi)
ghs=1/(1-xi(3))+(di*dj/(di+dj))*3*xi(2)/((1-xi(3))^2)+(((di*dj/(di+dj)))^2)*2*((xi(2))^2)/((1-xi(3))^3);
end