function ind =  ntoi(value,nrow,grid,n)

aux = min(value,min(grid(1,n)));
aux = max(aux,max(grid(1,1)));
step = (grid(1,n)-grid(1,1))/(n-1);

ind = round(((aux-grid(1,1)*ones(nrow,1))/step)+ones(nrow,1));
ind = ind';
end