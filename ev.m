%Value function

function ev_out = ev(cash,nrow,ncol,v,nco,prob,n,fy,eyp,grid,secondd,ret,reg_coef)
ev_out = 0.0;
ones= 1.0;

for ind1=1:n
  for ind2=1:n


     inc=(fy(1,ind1)*(eyp(1,ind2)+reg_coef*ret));
     inc2 = inc*ones+cash;
     inc2=min(inc2,grid(nco,1));
     inc2=max(inc2,grid(1,1));
     aux = splint(grid,v(:,1),secondd',nco,inc2,nrow,ncol);
     v1=(prob(1,ind1)*prob(1,ind2))'*aux;

     ev_out=ev_out+v1;
  end
    
end 
ev_out = ev_out';
end 

