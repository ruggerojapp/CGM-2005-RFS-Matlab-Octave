   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                                                                      %
   %       Paper: "Consumption and Portfolio Choice over the Life-Cycle", %
   %       Joao Cocco, Francisco Gomes and Pascal Maenhout,               %
   %       The Review of Financial Studies, 18 (2), 491-533, 2005.        %
   %                                                                      %
   %                        Replication by: Ruggero Jappelli              %
   %                                                                      %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc;
tic

%% Parameters
tb=20; tr=65; td=100;

nqp=3; nalfa=101; ncash=401; nc=1501;

delta=0.96; gamma= 10.0; 
infinity= -1e+10; rf=1.02;
sigma_r=0.157^2; exc= 0.04;
mu = 0; reg_coef=0.0;

ret_y=0; ret_fac=0.68212;

a=-2.170042+2.700381; b1=0.16818;
b2=-0.0323371/10; b3=0.0019704/100;

sigt_y=0.0738; sigp_y=0.01065;

sovrprob = zeros(1,td-tb);
delta2   = zeros(1,td-tb);

grid = ones(1,nqp); weig = ones(1,nqp);
gr   = ones(1,nqp); eyp  = ones(1,nqp); 
eyt  = ones(1,nqp); gret  = ones(1,nqp);
ones_nqp_1 = ones(1,nqp);

expeyp = zeros(1,nqp); f_y = zeros(tr-1,nqp);
gcash = zeros(ncash,1); ut = zeros(ncash,1);

aux3 = zeros(1,ncash); secd = zeros(1,ncash);

gc = zeros(nc,1); u = zeros(nc,1);

galfa = zeros(1,nalfa); c = zeros(ncash,2);
alfa = zeros(ncash,2);  v = zeros(ncash,2);

%% Quadrature 

 weig(1,1)= 0.1666666666666;
 weig(1,2)= 0.6666666666666;
 weig(1,3)= 0.1666666666666;
 grid(1,1)= -1.73205080756887;
 grid(1,2)=  0.0;
 grid(1,3)=  1.73205080756887;

 %% Conditional survival probabilities 
survprob(1,1) = 0.99845;
survprob(2,1) = 0.99839;
survprob(3,1) = 0.99833;
survprob(4,1) = 0.9983;
survprob(5,1) = 0.99827;
survprob(6,1) = 0.99826;
survprob(7,1) = 0.99824;
survprob(8,1) = 0.9982;
survprob(9,1) = 0.99813;
survprob(10,1) = 0.99804;
survprob(11,1) = 0.99795;
survprob(12,1) = 0.99785;
survprob(13,1) = 0.99776;
survprob(14,1) = 0.99766;
survprob(15,1) = 0.99755;
survprob(16,1) = 0.99743;
survprob(17,1) = 0.9973;
survprob(18,1) = 0.99718;
survprob(19,1) = 0.99707;
survprob(20,1) = 0.99696;
survprob(21,1) = 0.99685;
survprob(22,1) = 0.99672;
survprob(23,1) = 0.99656;
survprob(24,1) = 0.99635;
survprob(25,1) = 0.9961;
survprob(26,1) = 0.99579;
survprob(27,1) = 0.99543;
survprob(28,1) = 0.99504;
survprob(29,1) = 0.99463;
survprob(30,1) = 0.9942;
survprob(31,1) = 0.9937;
survprob(32,1) = 0.99311;
survprob(33,1) = 0.99245;
survprob(34,1) = 0.99172;
survprob(35,1) = 0.99091;
survprob(36,1) = 0.99005;
survprob(37,1) = 0.98911;
survprob(38,1) = 0.98803;
survprob(39,1) = 0.9868;
survprob(40,1) = 0.98545;
survprob(41,1) = 0.98409;
survprob(42,1) = 0.9827;
survprob(43,1) = 0.98123;
survprob(44,1) = 0.97961;
survprob(45,1) = 0.97786;
survprob(46,1) = 0.97603;
survprob(47,1) = 0.97414;
survprob(48,1) = 0.97207;
survprob(49,1) = 0.9697;
survprob(50,1) = 0.96699;
survprob(51,1) = 0.96393;
survprob(52,1) = 0.96055;
survprob(53,1) = 0.9569;
survprob(54,1) = 0.9531;
survprob(55,1) = 0.94921;
survprob(56,1) = 0.94508;
survprob(57,1) = 0.94057;
survprob(58,1) = 0.9357;
survprob(59,1) = 0.93031;
survprob(60,1) = 0.92424;
survprob(61,1) = 0.91717;
survprob(62,1) = 0.90922;
survprob(63,1) = 0.90089;
survprob(64,1) = 0.89282;
survprob(65,1) = 0.88503;
survprob(66,1) = 0.87622;
survprob(67,1) = 0.86576;
survprob(68,1) = 0.8544;
survprob(69,1) = 0.8423;
survprob(70,1) = 0.82942;
survprob(71,1) = 0.8154;
survprob(72,1) = 0.80002;
survprob(73,1) = 0.78404;
survprob(74,1) = 0.76842;
survprob(75,1) = 0.75382;
survprob(76,1) = 0.73996;
survprob(77,1) = 0.72464;
survprob(78,1) = 0.71057;
survprob(79,1) = 0.6961;
survprob(80,1) = 0.6809;
survprob(81,1) = 0;


survprob = survprob';
delta2 = delta*survprob;

%% Additional computations

tn = td-tb+1;

gr=grid*sigma_r^0.5;
eyp=grid*sigp_y^0.5;
eyt=grid*sigt_y^0.5;
mu = exc+rf;

expeyp = exp(eyp);

%% Construct grids

for ind1=1:nalfa
   galfa(1,ind1)=(ind1-1.0)/(nalfa-1.0);
end
gret = mu*ones_nqp_1+gr;
for ind1=1:ncash
   gcash(ind1,1)=4.0+(ind1-1.0)*1.0;
end
aux3(1,:) = gcash(:,1);
for ind1=1:nc
   gc(ind1,1)=0.0+(ind1-1.0)*0.25;
end


%% Labor Income 

for ind1=tb+1:tr
   avg = exp(a+b1*ind1+b2*ind1^2+b3*ind1^3);
   f_y(ind1-tb,:) = avg*exp(eyt(1,:));
end


ret_y= ret_fac*avg;

%% Terminal Period
 
ut = utility(gcash,ncash,gamma);

v(:,1)= ut(:,1);
c(:,1)= gcash(:,1);

%% Retirement Periods
u = utility(gc,nc,gamma);
tt=80;

disp("Simulating period")


for ind1= 1:35
t= tt-ind1+1;
disp(t)
secd = spline(aux3',v(:,1),ncash,gamma);


 for ind2= 1:ncash 

    if (t == tn-1) 
         lowc= c(ind2,1)/2.0;
         highc= c(ind2,1);
         if (gcash(ind2,1)>=50)
             highc= c(ind2,1)/1.5;
         end
    elseif t==tn-2
         lowc= c(ind2,1)/2.5;
         highc= c(ind2,1);
         if gcash(ind2,1)>=50
             highc= c(ind2,1)/1.2;
         end
    elseif (t<tn-2 && t>tn-5)
         lowc= c(ind2,1)/3.5;
         highc= c(ind2,1)+0.0;
         if (gcash(ind2,1)>=50)
            highc= c(ind2,1)/1.1;
         end
    else
         lowc= c(ind2,1)-10.0;
         highc= c(ind2,1)+10.0;
    end
 
    
    lowc2  = ntoi(lowc,1,gc',nc);
    highc2 = ntoi(highc,1,gc',nc);
   
    nc_r= highc2-lowc2+1;
    gc_r=zeros(nc_r,1);
    gc_r(:,1)= gc(lowc2:highc2,1);
    
    lowalfa2= 1.0;
    highalfa2= nalfa;
    
        if (gcash(ind2,1)>40 && t<tn-1)
          lowalfa   = alfa(ind2,1)-0.2;
          highalfa  = alfa(ind2,1)+0.2;
          lowalfa2  = ntoi(lowalfa,1,galfa,nalfa);
          highalfa2 = ntoi(highalfa,1,galfa,nalfa);
        end        
        
       nalfa_r= highalfa2-lowalfa2+1;
       galfa_r=zeros(1,nalfa_r);
       galfa_r(1,:) = galfa(1,lowalfa2:highalfa2);
       invest=zeros(1,nc_r);
       idvec = ones(1,nc_r);
       invest(1,:) = gcash(ind2,1)*idvec(:,1)-gc_r(:,1);
       
       u_r=zeros(1,nc_r);
       u_r(1,:) = u(lowc2:highc2,1);
       u2=zeros(1,nc_r);
       for ind4=1:nc_r
          if invest(1,ind4)<0
            u2(1,ind4) = infinity;
          else
            u2(1,ind4) = u_r(1,ind4);
          end
       end
               
       invest = max(invest,0);
       u3=zeros(nalfa_r,nc_r);
       for ind4=1:nalfa_r
          u3(ind4,:)=u2(1,:);
       end
       u3 = max(u3,infinity);
       v1=zeros(nalfa_r,nc_r);
       
       nw=zeros(nalfa_r,nc_r);
       nv=zeros(nalfa_r,nc_r);
       for ind5=1:nqp
          nw = fci(invest,nc_r,galfa_r,nalfa_r,gret(1,ind5),rf);
          nv = evr(nw,nc_r,nalfa_r,v(:,1),1,ncash,ret_y,aux3',secd);
          v1 = v1+nv*weig(1,ind5);
          
       end
       
         reslt = v1; 
       
       
        vv= zeros(nalfa_r,nc_r);
        vv = u3+delta2(1,t)*v1;
        
        
         
        vv = max(vv,infinity);
        auxv=zeros(1,nc_r*nalfa_r);
        
        %rj: FORTRAN stacks rows horizontally
        auxv = reshape(vv',nc_r*nalfa_r,1);
            
        
       
        
        [M,I] = max(auxv);
        v(ind2,2) = M;
        
        v(ind2,2);
        
        pt = I;
        aux2 = floor(real(pt-1)/real(nc_r));
        
        alfa(ind2,2) = galfa(1,aux2+lowalfa2);
        c(ind2,2) = gc(pt(1)-aux2*nc_r+lowc2-1,1);
        
       
 end

%Store Results  
  for ind5=1:ncash
    result(ind5,1,t) = alfa(ind5,2);
    result(ind5,2,t) = c(ind5,2);
    result(ind5,3,t) = v(ind5,2);
  end
  
  
v(:,1)=v(:,2);
c(:,1)=c(:,2);
alfa(:,1)=alfa(:,2); 
end



%% Other Periods
%Same story different value function -- EV instead of EVR (retirement)

for ind1= 1:tt-35
    t= 45-ind1+1;
    disp(t)
    secd(1,:) = spline(aux3',v(:,1),ncash,gamma);
   
    for ind2= 1:ncash
        if (t<tr-19 && t>tr-25) 
 	    lowc= c(ind2,1)-10.0;
 	    highc= c(ind2,1)+10.0;
        else
 	    lowc= c(ind2,1)-5.0;
 	    highc= c(ind2,1)+5.0;
    end

    lowc2 = ntoi(lowc,1,gc',nc);
    highc2 = ntoi(highc,1,gc',nc);

    nc_r= highc2-lowc2+1;
    gc_r=zeros(nc_r,1);
    gc_r(:,1)= gc(lowc2:highc2,1);
    lowalfa2= 1.0;
    highalfa2= nalfa;
        
    if (gcash(ind2,1)>40 && t<tn-1)
        lowalfa   = alfa(ind2,1)-0.2;
        highalfa  = alfa(ind2,1)+0.2;
        lowalfa2  = ntoi(lowalfa,1,galfa,nalfa);
        highalfa2 = ntoi(highalfa,1,galfa,nalfa);
    end   

    nalfa_r= highalfa2-lowalfa2+1;
    galfa_r=zeros(1,nalfa_r);
    galfa_r(1,:) = galfa(1,lowalfa2:highalfa2);
    invest=zeros(1,nc_r);
    idvec = ones(1,nc_r);
    invest(1,:) = gcash(ind2,1)*idvec(:,1)-gc_r(:,1);
       
    u_r=zeros(1,nc_r);
    u_r(1,:) = u(lowc2:highc2,1);
    u2=zeros(1,nc_r);
    for ind4=1:nc_r
        if invest(1,ind4)<0
            u2(1,ind4) = infinity;
        else
            u2(1,ind4) = u_r(1,ind4);
        end
    end
       
    invest = max(invest,0);
    u3=zeros(nalfa_r,nc_r);
       for ind4=1:nalfa_r
          u3(ind4,:)=u2(1,:);
       end
       
    u3 = max(u3,infinity);
    v1=zeros(nalfa_r,nc_r);
    nw=zeros(nalfa_r,nc_r);
    nv=zeros(nalfa_r,nc_r);
    
    for ind5=1:nqp
          nw = fci(invest,nc_r,galfa_r,nalfa_r,gret(1,ind5),rf);
          nv = ev(nw,nc_r,nalfa_r,v(:,1),ncash,weig,nqp,f_y(t,:),expeyp,aux3',secd,gret(1,ind5),reg_coef);
          v1 = v1+nv*weig(1,ind5);
    end

        reslt = v1;
        vv= zeros(nalfa_r,nc_r);
        vv = u3+delta2(1,t)*v1;
             
        vv = max(vv,infinity);
        auxv=zeros(1,nc_r*nalfa_r);
        
        %rj: FORTRAN stacks rows horizontally
        auxv = reshape(vv',nc_r*nalfa_r,1);
            
        [M,I] = max(auxv);
        v(ind2,2) = M;
        
        pt = I;
        aux2 = floor(real(pt-1)/real(nc_r));
        
        alfa(ind2,2) = galfa(1,aux2+lowalfa2);
        c(ind2,2) = gc(pt(1)-aux2*nc_r+lowc2-1,1);
        
    end  
  %Store Results
  for ind5=1:ncash
    result(ind5,1,t) = alfa(ind5,2);
    result(ind5,2,t) = c(ind5,2);
    result(ind5,3,t) = v(ind5,2);
  end

v(:,1)=v(:,2);
c(:,1)=c(:,2);
alfa(:,1)=alfa(:,2);
end      
disp("Done!")

%% Plot Results


figure
x = (1:1:300)';
c20 = result(1:300,2,1);
c35 = result(1:300,2,16);
c65 = result(1:300,2,46);
c85 = result(1:300,2,66);

plot(x, c20,x, c35,'--',x, c65,':', x, c85, "-")
title('Optimal Consumption')
xlabel('Cash-on-hand')
ylabel('Thousands of 1992 US dollars')
legend("Year 20", "Year 35", "Year 65", "Year 85")
legend('Orientation','horizontal', 'Location', 'southoutside')
%legend('boxoff')
saveas(gcf,"Consumption profile.pdf");



figure
x = (20:1:300)';
a30 = result(20:300,1,11);
a55 = result(20:300,1,36);
a75 = result(20:300,1,56);

%Note: the first α values are not shown in cgm (2005), they
%      serve to initialize the numerical optimization. 
%      See FORTRAN folder.

plot(x, a30,'--',x, a55,':', x, a75, "-")
title('Optimal Investment')
xlabel('Cash-on-hand')
ylabel('Portfolio share in stocks')
legend("Year 30", "Year 55", "Year 75")
legend('Orientation','horizontal', 'Location', 'southoutside')
%legend('boxoff')
saveas(gcf,"Investment profile.pdf");

