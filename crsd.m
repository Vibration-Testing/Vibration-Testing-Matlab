function [fout,Pxyout]=crsd(x,y,dt,n,ave)
%CRSD One sided cross spectrum density estimate of the data sequences X and Y (Sxy).
% [F,Pxy] = CRSD(X,Y,DT,N,AVE) estimates the Cross Spectrum
% Density using the direct method. 
% DT is the time between date samples.
% If DT is the time vector, DT is extracted as T(2)-T(1).
% N is the number of points to be used in the Fourier 
% Transform. The default for N is to zero pad the data 
% to the next power of 2.
% If X and Y are matrices, CRSD will find the Cross Spectrum 
% Density for each column and average the results
% unless AVE is set to 'noave'. N and AVE are optional.
% Either can be left out.
%
% CRSD(X,Y,DT,...) plots the Cross Spectrum Density if there
% are no output arguments.  Click in the region of interest
% to zoom in.  Each click will double the size of the plot.
% Double click to return to full scale.
%
% See page 38, R.K. Otnes and L. Enochson, Digital Time Series
% Analysis, J. Wiley and Sons, 1972.
%
% See also TFEST, ASD, TFPLOT

%	Copyright (c) 1994 by Joseph C. Slater
% Amplitude adjusted to match Vibration Testing... 3/31/2003
% Frequency Scale Fixed to work with a large number of 
%       short data sets. 12/27/98
% Normalized for magnitude 9/23/98
% Frequency Scale Fixed 9/23/98
% Fixed amplitude of zero frequency value 3/30/16 (filed by Admir Makas)

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
n=length(x);
ffty=fft(y,n)*dt;
fftx=fft(x,n)*dt;
Pxy=ffty.*conj(fftx)/(n*dt)*2;

Pxy(:,1) = Pxy(:,1)/2;


Pxy=Pxy(1:ceil(n/2),:);
if sy(2)~=1 & ~strcmp(ave,'noave')
 Pxy=mean(Pxy')';
end
nfreq=1/dt/2;
sfft=size(fftx);
lfft=sfft(1);
%f=(0:(nfreq/(lfft/2-1)):nfreq)';
fmax=(nfreq/(lfft/2))*(size(Pxy,1)-1);
%f=(0:(nfreq/(lfft/2)):nfreq-1)';
f=(0:(nfreq/(lfft/2)):fmax)';
if nargout==0
 semilogy(f,abs(Pxy))
 
 title('Cross Spectrum Density')
 xlabel('Frequency (Hz)')
 ylabel('Cross Spectrum Density')
 grid
 zoom on
 return
end

fout=f;
Pxyout=Pxy;
