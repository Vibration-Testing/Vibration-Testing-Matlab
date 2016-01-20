function f = expwin(n,ts)
%EXPWIN Returns or applies an N-point exponential window
% X=EXPWIN(N,TS) produces an N-point exponential window
%
% YWIN=EXPWIN(Y,TS) Returns Y with a exponential window
% applied to the vector Y or set of vectors represented 
% by the matrix Y.
%
% TS is the point in the window at which the exponential 
% reaches .05 (TS between 0 and 1). The default value 
% for TS is 0.75.
%
% See also BLACKWIN, BOXWIN, HAMMWIN, TRIWIN, and VONHANN.


%	Copyright (c) 1994 by Joseph C. Slater
% 4/6/2003 Normalization for ASD (accurate in limit)

if nargin==1
  ts=.75;
  disp('Using 75% window length for 5% window magnitude.')
end

if ts > 1 | ts < 0
 tstxt=num2str(ts);
 text=['Error: TS should be between 0 and 1. '...
       'TS was entered as ' tstxt '.'];
 error(text)
end

tc=ts/2.9957;

sn=size(n);

if sn==[1 1]
  v=(n-1)/n*(0:n-1)+(n-1)/n/2;
  f=exp(-v/tc/(n-1))';
  f=f/sqrt(f'*f)*sqrt(length(f));    
else
  %  expomesh=meshdom(expwin(sn(1),.75),1:sn(2))';  
  expomesh=meshgrid(expwin(sn(1),.75),1:sn(2))';
  f=n.*expomesh;
end

