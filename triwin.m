function f = triwin(n)
%TRIWIN Returns of applies Bartlett (triangular) window.
% X=TRIWIN(N) produces an N-point triangle window.
%
% YWIN=TRIWIN(Y) Returns Y with a triangle window
% applied to the vector Y or set of vectors 
% represented by the matrix Y.
%
% See also BLACKWIN, BOXWIN, EXPWIN, HAMMWIN, and 
% VONHANN.

%	Copyright (c) 1994 by Joseph C. Slater
% 4/6/2003 Normalization for ASD (accurate in limit)
  
sn=size(n);

if sn==[1 1]
  if n/2~=floor(n/2)
    v=((-n/2+.5):1:(n/2-.5))*2/n;
    f=1-abs(v');
    f=f/sqrt(f'*f)*sqrt(length(f));    
    %f=f*sqrt(3);
    %f=f*(sqrt(n/(n+2/n)));
  else
    v=(((-n/2)+.5):1:((n/2)-.5))/(n)*2;
    f=1-abs(v');
    %f=f*sqrt(3);
    f=f/sqrt(f'*f)*sqrt(length(f));
    %f=f*(sqrt(n/(n+2/n)));
  end
 else
  trimesh=meshdom(triwin(sn(1)),1:sn(2))';
  f=n.*trimesh;
end

