function [fr,ms]=soeig(m,k,n)
%  [fr,ms]=SOEIG(m,k,n) returns the first n natural 
%          frequencies (in Hz) and mode shapes of 
%          the second order system defined by m and k.

%  Copyright Joseph C. Slater, 1996
%  All rights reserved.
%  Added to Vibration Toolbox 9/23/98

l=max(size(m));
if nargin==2
  n=l;
end
if n>l
  disp(['Only ' num2str(l) ' exist.'])
  break
end
m=sparse(m);
k=sparse(k);
r=chol(m);
%whos;
kt=(r')\k/r;
kt=(kt+kt')/2;
[v,d]=eig(full(kt));
[d,i]=sort(sqrt(diag(d)/2/pi));
u=r\sparse(v);
u=u(:,i); 

fr=d(1:n);
ms=full(u(:,1:n));
