function y = splint(xa,ya,y2a,n,x,nrow,ncol)
klo = 1;
khi = n;

for indr=1:nrow
  for indc=1:ncol
       klo = 1;
       khi = n;
       while (khi-klo>1)
          k = round((khi+klo)/2);
          if (xa(k,1)>x(indr,indc))
            khi=k;
          else
            klo=k;
          end 
       end           

    h = xa(khi,1)-xa(klo,1);
    a = (xa(khi,1)-x(indr,indc))/h;
    b = (x(indr,indc)-xa(klo,1))/h;
       
    y(indr,indc) = a*ya(klo,1)+b*ya(khi,1)+((a^3-a)*y2a(klo,1)+(b^3-b)*y2a(khi,1))*(h^2)/6.0;
  end 
end 
