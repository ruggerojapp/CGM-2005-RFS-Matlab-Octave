%Value function during retirement

function ev_out = evr(cash,nrow,ncol,v,nro,nco,fy,grid,secondd)



ons = repmat(1,nrow,ncol);
inc=fy*ones+cash;

inc=min(inc,grid(nco,1));
inc=max(inc,grid(1,1));


prob_li = 0.0;
aux = splint(grid,v(:,1),secondd',nco,inc,nrow,ncol);
ev_out=aux';
end 
