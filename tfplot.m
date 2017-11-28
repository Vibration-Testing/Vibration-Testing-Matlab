function TFPLOT(F,Xfer,Fmin,Fmax,FLAG)
%TFPLOT Plots various transfer functions.
% TFPLOT(F,Xfer,Fmin,Fmax) plots the transfer function 
% between Fmin and Fmax. F is the frequency vector (in Hz),
% Xfer is the transfer function. Fmin and Fmax are the 
% minimum and maximum frequencies to be shown on the plot.
%
% TFPLOT(F,Xfer) plots the transfer function over the 
% entire frequency range of F.
%
% TFPLOT(F,Xfer,Fmin,Fmax,FLAG) plots the transfer 
% function in different forms depending on the value
% of FLAG.  Click in the region of interest
% to zoom in.  Each click will double the size of the plot.
% Double click to return to full scale.
%
%  FLAG(1)     Plot Type
%  -------     ---------
%    1 (def)   Magnitude and Phase versus F
%    2         Magnitude and Phase versus log10(F)
%    3         Bodelog  (Magnitude and Phase versus log10(w))
%    4         Real and Imaginary
%    5         Nyquist  (Real versus Imaginary)
%    6         Magnitude versus F 
%    7         Phase versus F
%    8         Real versus F
%    9         Imaginary versus F
%   10         Magnitude versus log10(F) 
%   11         Phase versus log10(F)
%   12         Real versus log10(F)
%   13         Imaginary versus log10(F)
%   14         Magnitude versus log10(w) 
%   15         Phase versus log10(w)       
% 
%
%  TFPLOT(F,Xfer,FLAG) plots the transfer function 
%  over the entire frequency range of F in the 
%  form defined by FLAG.
%
%  Example (copy and paste into Matlab command line):
%  f=(0:.01:100)';
%  w=f*2*pi;
%  k=1e5;m=1;c=1;
%  tf=1./(m*(w*j).^2+c*j*w+k);
%  figure(1);tfplot(f,tf)
%  figure(2);tfplot(f,tf,5)
%
%    See also TFEST, ASD, and CRSD.

% Copyright J. Slater, Dec 17, 1994
% Updated April 27, 1995

%Copy protection, expires Jan 1, 2001
lenF=length(F);

if lenF==1;
 F=(0:length(Xfer)-1)'*F;
end
if nargin==2
  Fmin=min(F);
  Fmax=max(F);
  FLAG=1;
 elseif nargin==3
  FLAG=Fmin;
  if FLAG~=2
    Fmin=min(F);
    Fmax=max(F);
   else
    Wmin=min(F)*2*pi;
    Wmax=max(F)*2*pi;
  end
 elseif nargin==4
  Fmax=min([max(F) Fmax]);
  Fmin=max([min(F) Fmin]);
  FLAG=1;
 else
  if (FLAG~=3 & FLAG~=14 & FLAG~=15)
    Fmax=min([max(F) Fmax]);
    Fmin=max([min(F) Fmin]);  
   else
    Wmax=min([max(F)*2*pi Fmax*2*pi]);
    Wmin=max([min(F)*2*pi Fmin*2*pi]);  
  end
end

if Fmin>Fmax
  disp('Fmin must be less than Fmax.')
  return
end


inlow=floor(length(F)*(Fmin-min(F))/(max(F)-min(F)))+1;

inhigh=ceil(length(F)*(Fmax-min(F))/(max(F)-min(F)))+1;
if inlow<1,inlow=1;end
if inhigh>lenF,inhigh=lenF;end
Xfer=Xfer(inlow:inhigh,:);
F=F(inlow:inhigh,:);
mag=20*log10(abs(Xfer));
mag;
minmag=min(min(mag));
maxmag=max(max(mag));
phase=unwrap(angle(Xfer))*180/pi;
phmin_max=[floor(min(min(phase))/45)*45 ceil(max(max(phase))/45)*45];
minreal=min(min(real(Xfer)));
maxreal=max(max(real(Xfer)));
minimag=min(min(imag(Xfer)));
maximag=max(max(imag(Xfer)));


if FLAG==1
  subplot(2,1,1)
  plot(F,mag)
  xlabel('Frequency (Hz)')
  ylabel('Mag (dB)')
  grid on
%  Fmin,Fmax,min(mag),max(mag)
  axis([Fmin Fmax minmag maxmag])
  zoom on
  subplot(2,1,2)
  plot(F,phase)
  xlabel('Frequency (Hz)')
  ylabel('Phase (deg)')
  grid on
  axis([Fmin Fmax  phmin_max(1) phmin_max(2)])
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
  zoom on
 elseif FLAG==2
  subplot(2,1,1)
  semilogx(F,mag)
  xlabel('Frequency (Hz)')
  ylabel('Mag (dB)')
  grid on
%  Fmin,Fmax,min(mag),max(mag)
  axis([Fmin Fmax minmag maxmag])
  zoom on
  subplot(2,1,2)
  semilogx(F,phase)
  xlabel('Frequency (Hz)')
  ylabel('Phase (deg)')
  grid on
  axis([Fmin Fmax  phmin_max(1) phmin_max(2)])
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
  zoom on
 elseif FLAG==3
  subplot(2,1,1)
  mag=20*log10(abs(Xfer));
  semilogx(F*2*pi,mag)
  xlabel('Frequency (Rad/s)')
  ylabel('Mag (dB)')
  grid on
  axis([Wmin Wmax minmag maxmag])
  zoom on
  subplot(2,1,2)
  semilogx(F*2*pi,phase)
  xlabel('Frequency (Rad/s)')
  ylabel('Phase (deg)')
  grid on
  axis([Wmin Wmax  phmin_max(1) phmin_max(2)])
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
  zoom on
 elseif FLAG==4
  subplot(2,1,1)
  plot(F,real(Xfer))
  xlabel('Frequency (Hz)')
  ylabel('Real')
  grid on
  axis([Fmin Fmax minreal maxreal])
  zoom on
  subplot(2,1,2)
  plot(F,imag(Xfer))
  xlabel('Frequency (Hz)')
  ylabel('Imaginary')
  grid on
  axis([Fmin Fmax minimag maximag])
  zoom on
 elseif FLAG==5
  subplot(1,1,1)
  imax=round(length(F)*Fmax/max(F));
  imin=round(length(F)*Fmin/max(F))+1;
  plot(real(Xfer(imin:imax)),imag(Xfer(imin:imax)))
  xlabel('Real')
  ylabel('Imaginary')
  grid on
  zoom on
 elseif FLAG==6
  subplot(1,1,1)
  mag=20*log10(abs(Xfer));
  plot(F,mag)
  xlabel('Frequency (Hz)')
  ylabel('Mag (dB)')
  grid on
  axis([Fmin Fmax minmag maxmag])
  zoom on
 elseif FLAG==7
  subplot(1,1,1)
  plot(F,phase)
  xlabel('Frequency (Hz)')
  ylabel('Phase (deg)')
  grid on
  phmin_max=[floor(min(phase)/45)*45 ceil(max(phase)/45)*45];
  axis([Fmin Fmax  phmin_max(1) phmin_max(2)])
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
  zoom on
 elseif FLAG==8
  subplot(1,1,1)
  plot(F,real(Xfer))
  xlabel('Frequency (Hz)')
  ylabel('Real')
  grid on
  axis([Fmin Fmax minreal maxreal])
  zoom on
 elseif FLAG==9
  subplot(1,1,1)
  plot(F,imag(Xfer))
  xlabel('Frequency (Hz)')
  ylabel('Imaginary')
  grid on
  axis([Fmin Fmax minimag maximag])
  zoom on
 elseif FLAG==10
  subplot(1,1,1)
  mag=20*log10(abs(Xfer));
  semilogx(F,mag)
  xlabel('Frequency (Hz)')
  ylabel('Mag (dB)')
  grid on
  axis([Fmin Fmax minmag maxmag])
  zoom on
 elseif FLAG==11
  subplot(1,1,1)
  semilogx(F,phase)
  xlabel('Frequency (Hz)')
  ylabel('Phase (deg)')
  grid on
  phmin_max=[floor(min(phase)/45)*45 ceil(max(phase)/45)*45];
  axis([Fmin Fmax  phmin_max(1) phmin_max(2)])
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
  zoom on
 elseif FLAG==12
  subplot(1,1,1)
  semilogx(F,real(Xfer))
  xlabel('Frequency (Hz)')
  ylabel('Real')
  grid on
  axis([Fmin Fmax minreal maxreal])
  zoom on
 elseif FLAG==13
  subplot(1,1,1)
  semilogx(F,imag(Xfer))
  xlabel('Frequency (Hz)')
  ylabel('Imaginary')
  grid on
  axis([Fmin Fmax minimag maximag])
  zoom on
 elseif FLAG==14
  subplot(1,1,1)
  mag=20*log10(abs(Xfer));
  semilogx(F*2*pi,mag)
  xlabel('Frequency (Rad/s)')
  ylabel('Mag (dB)')
  grid on
  axis([Wmin Wmax minmag maxmag])
  zoom on
 elseif FLAG==15
  subplot(1,1,1)
  semilogx(F*2*pi,phase)
  xlabel('Frequency (Rad/s)')
  ylabel('Phase (deg)')
  grid on
  axis([Wmin Wmax  phmin_max(1) phmin_max(2)])
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
  zoom on
 else
  subplot(2,1,1)
  mag=20*log10(abs(Xfer));
  plot(F,mag)
  xlabel('Frequency (Hz)')
  ylabel('Mag (dB)')
  grid on
  axis([Fmin Fmax minmag maxmag])
  zoom on
  subplot(2,1,2)
  plot(F,phase)
  xlabel('Frequency (Hz)')
  ylabel('Phase (deg)')
  grid on
  phmin_max=[floor(min(phase)/45)*45 ceil(max(phase)/45)*45];
  axis([Fmin Fmax phmin_max(1) phmin_max(2)])
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
  zoom on
end
