function f = vonhann(n)
%VONHANN Returns or applies a von Hann window.
% X=VONHANN(N) produces an N-point von Hann window
%
% YWIN=VONHANN(Y) Returns Y with a von Hann window
% applied to the vector Y or set of vectors 
% represented by the matrix Y.
%
% The von Hann window is more commonly known as 
% the "hanning" window and derives its name from 
% Julius von Hann.
%
% See also BLACKWIN, BOXWIN, EXPWIN, HAMMWIN, and TRIWIN.

% Copyright (c) 1994 by Joseph C. Slater
% Updated July 20, 1998 to use meshgrid
sn=size(n);

if sn==[1 1]
  f=(sin(pi*((0:n-1)+.5)/(n)).^2)'*sqrt(8/3);
 else
  hannmesh=meshgrid(hannwin(sn(1)),1:sn(2))';
  f=n.*hannmesh;
end
