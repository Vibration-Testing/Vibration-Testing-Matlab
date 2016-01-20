function [fout,Pyyout]=asd(y,dt,n,ave)
%ASD Onesided Auto or power spectrum density estimate of the data sequence y.
% [F,Pyy] = ASD(Y,DT,N,AVE) estimates the Power/Auto Spectral
% Density using the direct method. N is the number of 
% points to be used in the Fourier Transform. The default 
% is not to zero pad to the next power of 2. If DT is the time 
% vector, DT is extracted as T(2)-T(1). If Y is a matrix, 
% ASD will find the Auto Spectral Density for each column 
% and average the results unless AVE is set to 'noave'. 
% N and AVE are optional. Either can be left out.
%
% Zero padding is turned off harshly now. 
%
% ASD(Y,DT,N,AVE) plots the Auto Spectral Density if there
% are no output arguments. Click in the region of interest
% to zoom in.  Each click will double the size of the plot.
% Double click to return to full scale.
%
% See page 38, R.K. Otnes and L. Enochson, Digital Time Series
% Analysis, J. Wiley and Sons, 1972.
%
% See also TFEST, CRSD, TFPLOT

% Copyright (c) 1994 by Joseph C. Slater
% Amplitude adjusted to match Vibration Testing...3/31/2003
% Frequency Scale Fixed to work with a large number of 
%       short data sets. 12/27/98
% Frequency Scale Fixed 7/21/98
% Normalized for magnitude 7/21/98

n=length(y);
sy=size(y);
if nargin==2
  %n=2^nextpow2(sy(1));
  
  ave='yes';
 elseif nargin==3
  if strcmp(n,'noave')
   ave=n;
   n=length(y);
  else
   ave='yes';
  end
end

if length(dt)>1
 dt=dt(2)-dt(1);
end

if isempty(n)
  n=2^nextpow2(sy(1));
end
%size(y)
%n
%pause
ffty=fft(y,n)*dt;
%size(ffty)
%for ii=1:size(y,2)
%    ffty2(:,ii)=fft(y(:,ii),n);
%end
%size(ffty2)
%pause



Pyy=real(ffty.*conj(ffty))/(n*dt)*2;

Pyy=Pyy(1:ceil(n/2),:);
Pyy(1)=Pyy(1)/2;
%Pyy=Pyy((n:-1:n/2+1)+1,:);
if sy(2)~=1 & ~strcmp(ave,'noave')
 Pyy=mean(Pyy')';
end
nfreq=1/dt/2;
sfft=size(ffty);
lfft=sfft(1);
size(Pyy);
fmax=(nfreq/(lfft/2))*(size(Pyy,1)-1);
%f=(0:(nfreq/(lfft/2)):nfreq-1)';
f=(0:(nfreq/(lfft/2)):fmax)';
size(f);
if nargout==0
 plot(f,10*log10(Pyy))%see appendix
 
 title('Power Spectral Density of F(t) (dB)')
 xlabel('Frequency (Hz)')
 ylabel('Power Spectral Density')
 grid
 zoom on
 return
end

fout=f;
Pyyout=Pyy;
