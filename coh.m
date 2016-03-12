function [fout,cohout]=coh(x,y,dt,n,ave)
%COH Coherance.
% [F,COH]=COH(X,Y,DT,N,AVE) determines the coherance between 
% signals X and Y using CRSD and ASD.
% The coherance should be near 1 near the poles and may be
% near zero near the zeros without causing undue concern.
% DT is the time between samples. N is the number of 
% points to be used in the Fourier Transform.
% The default for N is to zero pad the data to the next 
% power of 2.
% If DT is the time vector, DT is extracted as T(2)-T(1).
% If X and Y are matrices, averaging will be performed on the
% power spectrums before the coherance is calculated
% unless AVE is set to 'noave'. N and AVE are optional.
% Either can be left out.
%
% COH(X,Y,DT,N) plots the Coherance if there are no ouput 
% arguments.  Click in the region of interest  to zoom in. 
% Each click will double the size of the plot. Double click 
% to return to full scale.
%
% See also TFEST, ASD, CRSD, and TFPLOT.

%	Copyright (c) 1994 by Joseph C. Slater
sy=size(y);
if nargin==3
  n=2^nextpow2(sy(1));
  ave='yes';
 elseif nargin==4
  if strcmp(n,'noave')
   ave=n;
   n=2^nextpow2(sy(1));
  else
   ave='yes';
  end
end

if isempty(n)
  n=2^nextpow2(sy(1));
end

if length(dt)>1
 dt=dt(2)-dt(1);
end

[f,Pxx]=crsd(x,x,dt,n,ave);
[f,Pxy]=crsd(x,y,dt,n,ave);
[f,Pyy]=crsd(y,y,dt,n,ave);
coh=abs(Pxy).^2./(Pxx.*Pyy);
if nargout==0
 plot(f,coh)
 %logo
 title('Coherance')
 xlabel('Frequency (Hz)')
 ylabel('Coherance')
 grid
 zoom on
 return
end

cohout=coh;
fout=f;
