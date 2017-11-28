function [a,p]=household(b,tol)
% HOUSEHOLD  Generates Householder matrix 
%    [A,P]=HOUSEHOLD(B) returns A and P such that
%    P'*B*P = A

% Copyright Joseph C. Slater, 5/2000
l=length(b)-2;
ii=eye(l+2);
pp=ii;
binit=b;
if nargin==1
	tol=eps;
end

if length(b)<1
	disp('This matrix is already tridiagonal.')
else

for i=1:l
	v(1,1:i)=zeros(1,i);
	alpha=norm(b(i,i+1:l+2))
	sn=sign(b(i,i+1));
	v(1,i+1)=sqrt((1+b(i,i+1)*sn/alpha)/2);
		for j=i+2:l+2
			v(1,j)=sn*b(i,j)/2/alpha/v(1,i+1);
		end
		v
	P=ii-2*v'*v
	pp=pp*P;
	b=pp'*binit*pp
end
a=b;
%a=pp'*binit*pp;
p=pp;

j=1:l;
er=max(abs(diag(a(j,j+2)))./abs(diag(a(j,j+1))))
if er>tol
	disp('Recursion')
	ainit=a;
	[a,pp]=household(a,tol)
	p=p*pp;	
	a=p'*binit*p
%	disp('Recursion')
	
end
end	
