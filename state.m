function [a,b,c,d]=state(m,d,k,bs)
% STATE [a,b,c,d]=STATE(M,C,K,Bs) returns the state matrix
% [xdot ]   [   0        I  ][ x  ]    [   0   ]
% [     ] = [               ][    ] +  [       ] [u]
% [xddot]   [-M^-1*K -M^-1*C][xdot]    [M^-1*Bs]
%
% from the originating equation
% Mxddot+Cxdot+x=Bs u
% 
% An identicty c matrix, and a zero d matrix are created.
  
% Copyright Joseph C. Slater, 2003
if nargin==2
  k=d;d=0*k;
end
l=length(diag(m));
a=[zeros(l) eye(l);-m\k -m\d];
len=size(m,1);
if nargin==3
  Bs=eye(len);
end

b=[zeros(len,size(Bs,2)); M\Bs];
c=eye(2*len);
d=zeros(size(C,1),size(B,2));
