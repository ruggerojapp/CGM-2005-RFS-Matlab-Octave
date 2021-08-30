function y2 = spline(x,y,n,gam)

yp1 = x(1,1)^(-gam);
y2(1,1)=-0.5;
u(1,1)=(3.0/(x(2,1)-x(1,1)))*((y(2,1)-y(1,1))/(x(2,1)-x(1,1))-yp1);

    for i=2:n-1
   sig = (x(i,1)-x(i-1,1))/(x(i+1,1)-x(i-1,1));
   p = sig*y2(i-1,1)+2.0;
   y2(i,1) = (sig-1.0)/p;
   u(i,1) = (6.0*((y(i+1,1)-y(i,1))/(x(i+1,1)-x(i,1))-(y(i,1)-y(i-1,1))/(x(i,1)-x(i-1,1)))/(x(i+1,1)-x(i-1,1))-sig*u(i-1,1))/p;
    end

y2(n,1) =0.0;

    for k=n-1:-1:1
   y2(k,1) = y2(k,1)*y2(k+1,1)+u(k,1);
    end
y2 = y2';
end 




















% function y2 = spline(x,y,n,gam)
% 
% yp1 = x(1,1)^(-gam);
% y2(1,1)=-0.5;
% u(1,1)=(3.0/(x(1,2)-x(1,1)))*((y(2,1)-y(1,1))/(x(1,2)-x(1,1))-yp1);
% for i=2:n-1
%    sig = (x(1,i)-x(1,i-1))/(x(1,i+1)-x(1,i-1));
%    p = sig*y2(i-1,1)+2.0;
%    y2(i,1) = (sig-1.0)/p;
%    u(i,1) = (6.0*((y(i+1,1)-y(i,1))/(x(1,i+1)-x(1,i))-(y(i,1)-y(i-1,1))/(x(1,i)-x(1,i-1)))/(x(1,i+1)-x(1,i-1))-sig*u(i-1,1))/p;
% end 
% y2(n,1) =0.0;
% for k=n-1:-1,1;
%    y2(k,1) = y2(k,1)*y2(k+1,1)+u(k,1);
% end 
% y2 = y2';
% end 


% implicit none
% INTEGER, INTENT(IN) :: n
% REAL, INTENT(IN) :: gam
% REAL, INTENT(IN) :: x(n,1), y(n,1)
% REAL, INTENT(OUT) :: y2(n,1)
% INTEGER :: i, k
% !REAL :: qn, un
% REAL :: yp1
% REAL :: p, sig
% REAL, ALLOCATABLE :: u(:,:)
% ALLOCATE(u(n,1))
% 
% 
