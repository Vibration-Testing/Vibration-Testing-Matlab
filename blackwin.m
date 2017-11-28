function f = blackwin(n)
%BLACKWIN Returns or applies a Blackman window.
% X=BLACKWIN(N) produces an N-point Blackman window
%
% YWIN=BLACKWIN(Y) Returns Y with a Blackman window
% applied to the vector Y or set of vectors 
% represented by the matrix Y.
%
% See also BOXWIN, EXPWIN, HAMMWIN, TRIWIN, and VONHANN.

% Copyright (c) 1994 by Joseph C. Slater

sn=size(n);

if sn==[1 1]
  f=(0.42-0.5*cos(2*pi*((0:n-1)+.5)/(n))+0.08*cos(4*pi*((0:n-1)+.5)/n))'*sqrt(5000/1523);
 else
  blackmesh=meshgrid(blackwin(sn(1)),1:sn(2))';
  f=n.*blackmesh;
end

