function f = hammwin(n)
%HAMMWIN Returns or applies a Hamming window.
% X=HAMMWIN(N) produces an N-point Hamming window
%
% YWIN=HAMMWIN(Y) Returns Y with a Hamming window
% applied to the vector Y or set of vectors 
% represented by the matrix Y.
%
% See also BLACKWIN, BOXWIN, EXPWIN, TRIWIN, and VONHANN.

% Copyright (c) 1994 by Joseph C. Slater

sn=size(n);

if sn==[1 1]
  f=(0.54-0.46*cos(2*pi*((0:n-1)+.5)/(n))')*sqrt(5000/1987);
 else
  hammmesh=meshgrid(hammwin(sn(1)),1:sn(2))';
  f=n.*hammmesh;
end
