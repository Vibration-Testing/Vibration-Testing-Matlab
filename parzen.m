function f = parzen(n)
%PARZWIN Returns or applies a Parzen window.
% X=PARZEN(N) produces an N-point Parzen window
%
% YWIN=PARZEN(Y) Returns Y with a Parzen window
% applied to the vector Y or set of vectors 
% represented by the matrix Y.
%
% See also BLACKWIN, BOXWIN, EXPWIN, TRIWIN, HAMMWIN, and VONHANN.
% Doesn't work for odd number of points.

% Copyright (c) 2003 by Joseph C. Slater
% 4/6/2003 Normalization for ASD (accurate in limit)
sn=size(n);

if sn==[1 1]
  normt=(1/(n-1):2/(n-1):1)';
  m=floor(n/4);
  normt1=normt(1:m);
  normt2=normt(m+1:floor(n/2));
  f1=1-6*normt1.^2+6*normt1.^3;
  f2=2*(1-normt2).^3;
  f=[f1; f2];
  f=[f(length(f):-1:1,1) ;f];
  f=f/sqrt(f'*f)*sqrt(length(f));    
 else
  parzmesh=meshgrid(parzwin(sn(1)),1:sn(2))';
  f=n.*parzmesh;
end
