function [tout,crcorout]=crcor(x,y,dt,type,ave)
%CRCOR Cross correlation.
% [Tau,COR]=CRCOR(X,Y,DT,TYPE,AVE) returns the Cross Correlation 
% between signals X and Y.
% [Tau,COR]=CRCOR(X,X,DT,TYPE,AVE) returns the Auto Correlation 
% of the signal X.
% DT is the time between samples.
% If DT is the time vector, DT is extracted as T(2)-T(1).
% TYPE is the type of correlation. TYPE = 1 causes CRCOR
% to return the linear correlation function. TYPE = 2
% causes CRCOR to return the circular correlation function.
% The default value is 1.
% If X and Y are matrices, averaging will be performed on the
% Correlations unless AVE is set to 'noave'. TYPE and AVE are 
% optional. Either can be left out.
%
% COH(X,Y,DT,N,AVE) plots the Correlation if there are no ouput 
% arguments. Click in the region of interest to zoom in. 
% Each click will double the size of the plot. Double click 
% to return to full scale.
%
% See also TFEST, ASD, COH, CRSD, and TFPLOT.

%	Copyright (c) 1994 by Joseph C. Slater
sy=size(y);
sy=size(y);
if nargin==3
  type=1;
  ave='yes';
 elseif nargin==4
  if strcmp(type,'noave')
   ave=n;
   type=1;
  else
   ave='yes';
  end
end

if isempty(type)
  type=1;
end

sx=size(x);
nc=sx(2);

if type==1
  n=sx(1)*2;
 else
  n=sx(1);
end

if length(dt)~=1
 dt=dt(2)-dt(1);
end

tmax=dt*(length(x)-1);
t=(-tmax:(2*tmax/(n-1)):tmax)'-(tmax/(n-1));

X=fft(x,n);
Y=fft(y,n);
pxy=real(ifft(conj(X).*Y));

crcr=fftshift(real(pxy));
crcr=crcr(1:length(crcr),:);

if nc~=1 & ~strcmp(ave,'noave')
 crcr=mean(crcr')';
end

if nargout==0
 plot(t,crcr)
 %logo
 if type==1
   text1='Linear ';
  else
   text1='Circular ';
 end
 if x==y
   text2='Auto ';
  else
   text2='Cross ';
 end
 text3=[text1 text2 'Correlation'];
 title(text3)
 xlabel('Time')
 ylabel(text3)
 grid
 zoom on
 return
end

crcorout=crcr;
tout=t;
