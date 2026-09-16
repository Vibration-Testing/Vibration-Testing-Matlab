function [fr,ms]=soeig(m,k,n)
%  [fr,ms]=SOEIG(m,k,n) returns the first n natural
%          frequencies (in Hz, not rad/s) and mass-normalized mode
%          shapes of the second order system defined by mass matrix m
%          and stiffness matrix k.
%
%  m and k must use consistent units (e.g., kg and N/m). fr(i) is the
%  i-th natural frequency in Hz; the corresponding angular natural
%  frequency in rad/s is 2*pi*fr(i).

%  Copyright Joseph C. Slater, 1996
%  All rights reserved.
%  Added to Vibration Toolbox 9/23/98

l=max(size(m));
if nargin==2
  n=l;
end
if n>l
  disp(['Only ' num2str(l) ' exist.'])
  n=l;
end
m=sparse(m);
k=sparse(k);
r=chol(m);
%whos;
kt=(r')\k/r;
kt=(kt+kt')/2;
[v,d]=eig(full(kt));
% d holds squared angular natural frequencies (rad/s)^2; convert to
% natural frequency in Hz via sqrt(d) (rad/s) / (2*pi).
[d,i]=sort(sqrt(diag(d))/2/pi);
u=r\sparse(v);
u=u(:,i); 

fr=d(1:n);
ms=full(u(:,1:n));
