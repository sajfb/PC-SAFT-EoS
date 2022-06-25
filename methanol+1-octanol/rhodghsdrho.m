function rhodghsdrho=rhodghsdrho(di,dj,xi)
rhodghsdrho=xi(3)/((1-xi(3))^2)+(di*dj/(di+dj))*(3*xi(2)/((1-xi(3))^2)+6*xi(2)*xi(3)/((1-xi(3))^3))+(((di*dj/(di+dj)))^2)*(4*((xi(2))^2)/((1-xi(3))^3)+6*((xi(2))^2)*xi(3)/((1-xi(3))^4));
end