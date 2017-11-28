function f = boxwin(n)
%BOXWIN Returns or applied a Boxcar window
% X=BOXWIN(N) produces an N-point Boxcar window
%
% YWIN=BOXWIN(Y) Returns Y with a Boxcar window
% applied to the vector Y or set of vectors 
% represented by the matrix Y.
%
% See also BLACKWIN, EXPWIN, HAMMWIN, TRIWIN, and VONHANN.

% Copyright (c) 1994 by Joseph C. Slater

sn=size(n);

if sn==[1 1]
  f=ones(n,1);
 else
  f=n;
end
