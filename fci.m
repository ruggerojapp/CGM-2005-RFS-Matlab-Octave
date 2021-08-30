%Portfolio returns

function capinc = fci(sav,nrow,galfa,n,ret,rf)
rp=ones(1,n);
for ind1=1:nrow
   rp = ret*galfa+rf*(ones(1,n)-galfa);
   for ind2=1:n
      capinc(ind1,ind2) = sav(1,ind1)*rp(1,ind2);
   end
end

end 