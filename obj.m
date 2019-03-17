%Implements Rosenbrock's function f
%If requested by calling function it calculated gradient g and hessian h
function [f,g,h]=obj(x);
n=max(size(x));
factor=1;
tot=(1.0 - x(1))^2;
for i= 2:n 
  factor= factor;
  tot= tot + factor*(x(i)-x(i-1)^2 )^2;
end;
f=tot;

if nargout > 1 
  g=zeros(n,1);
  factor=1;
  g(1)=2.0*(x(1)-1.0)-4.0*x(1)*(x(2)-x(1)^2 );
  h=zeros(n,n);
  h(1,1)=2+12*x(1)^2-4*x(2);
  h(1,2)=-4*x(1);
  h(2,1)=-4*x(1);
  for i= 2:n-1 
    factor= factor;
    g(i)=2.0*(x(i)-x(i-1)^2 )*factor -4.0*x(i)*(x(i+1)-x(i)^2 ...
						   )*factor;
    h(i,i)=2+12*x(i)^2-4*x(i+1);
    h(i,i+1)=-4*x(i);
    h(i+1,i)=-4*x(i);
  end;
  g(n)= 2.0*(x(n)-x(n-1)^2 )*factor;
  h(n,n)=2;
  df=g;




end;
