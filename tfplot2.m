function TFPLOT(F,Xfer,Fmin,Fmax,FLAG,junk)
%TFPLOT Plots various transfer functions.
% TFPLOT(F,Xfer,Fmin,Fmax) plots the transfer function 
% between Fmin and Fmax. F if the frequency vector (Hz),
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
%    1 (def)   Bodelin  (Magnitude and Phase versus F)
%    2         Bodelog  (Magnitude and Phase versus log10(W))
%    3         Real and Imaginary
%    4         Nyquist  (Real versus Imaginary)
%    5         Magnitude Only
%    6         Phase Only
%    7         Real Only
%    8         Imaginary Only
%
%  FLAG(1)=2 plots versus the log base 10 of the frequency in
%  Rad/s.
%
%  TFPLOT(F,Xfer,FLAG) plots the transfer function 
%  over the entire frequency range of F in the 
%  form defined by FLAG.
%
%    See also TFEST, ASD, and CRSD.

% Copyright J. Slater, Dec 17, 1994
if length(F)==1;
 F=(0:length(Xfer)-1)'*F;
end
if nargin==2
  Fmin=min(F);
  Fmax=max(F);
  FLAG=1;
 elseif nargin==3
  FLAG=Fmin;
  Fmin=min(F);
  Fmax=max(F);
 elseif nargin==4
  Fmax=min([max(F) Fmax]);
  Fmin=max([min(F) Fmin]);
  FLAG=1;
end
if nargin~=6
 global tfinfo
 gcf
 tfinfo(gcf,2:5)=[Fmin Fmax length(Xfer) F(2)-F(1)];
 tfinfo(gcf,1:length(Xfer)+5)=[tfinfo(gcf,1:5) Xfer'];
 size(tfinfo)
 tfinfo(:,1:8)
% tfinfo: Each row is for a new figure (new call to function)
% 1st column is for plot type
% 2nd column if for Fmin
% 3rd column is for Fmax
% 4th column is number of points in transfer function
% 5th column is delta F
% The rest are the data points of the transfer function
 tfinfo(gcf,1)=uicontrol(gcf,'style','popup','string',...
 'Bode|Imaginary|Magnitude|Nyquist|Phase|Real|Real and Imaginary',...
 'units','normalized',...
 'position',[.81 .95 .25 .05],...
 'callback',[...
 'global tfinfo ,',...
 'cb_col=[2,8,5,4,6,7,3];',...
 'tfinfo(gcf,1),',...
 'gcf,',...
 'tfinfo(gcf,4),'...
 'disp(''no. points''),',...
 'cb_col(get(tfinfo(gcf,1),''value'')),'...
 'tfplot2(tfinfo(gcf,5),tfinfo(gcf,6:5+tfinfo(gcf,4))'',tfinfo(gcf,2),tfinfo(gcf,3),cb_col(get(tfinfo(gcf,1),''value'')),1)']);
 format short
end
% 'cb_col=[''2'',''8'',''5'',''4'',''6'',''7'',''3''];',...

   
Fmax=min([max(F) Fmax]);
Fmin=max([min(F) Fmin]);

if FLAG==1
  subplot(2,1,1)
  mag=20*log10(abs(Xfer));
  plot(F,mag)
  xlabel('Frequency (Hz)')
  ylabel('Mag (dB)')
  grid on
  axis([Fmin Fmax min(mag) max(mag)])
  zoom on
  subplot(2,1,2)
  phase=unwrap(angle(Xfer))*180/pi;
  plot(F,phase)
  xlabel('Frequency (Hz)')
  ylabel('Phase (deg)')
  grid on
  phmin_max=[floor(min(phase)/45)*45 ceil(max(phase)/45)*45];
  axis([Fmin Fmax  phmin_max(1) phmin_max(2)])
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
  zoom on
 elseif FLAG==2
  subplot(2,1,1)
  mag=20*log10(abs(Xfer));
  semilogx(F*2*pi,mag)
  xlabel('Frequency (Rad/s)')
  ylabel('Mag (dB)')
  grid on
  Wmin=min(F)*2*pi;
  Wmax=max(F)*2*pi;
  axis([Wmin Wmax min(mag) max(mag)])
  zoom on
  subplot(2,1,2)
  phase=unwrap(angle(Xfer))*180/pi;
  semilogx(F*2*pi,phase)
  xlabel('Frequency (Rad/s)')
  ylabel('Phase (deg)')
  grid on
  phmin_max=[floor(min(phase)/45)*45 ceil(max(phase)/45)*45];
  axis([Fmin Fmax  phmin_max(1) phmin_max(2)])
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
  zoom on
 elseif FLAG==3
  subplot(2,1,1)
  plot(F,real(Xfer))
  xlabel('Frequency (Hz)')
  ylabel('Real')
  grid on
  axis([Fmin Fmax min(real(Xfer)) max(real(Xfer))])
  zoom on
  subplot(2,1,2)
  plot(F,imag(Xfer))
  xlabel('Frequency (Hz)')
  ylabel('Imaginary')
  grid on
  axis([Fmin Fmax min(imag(Xfer)) max(imag(Xfer))])
  zoom on
 elseif FLAG==4
  subplot(1,1,1)
  imax=round(length(F)*Fmax/max(F));
  imin=round(length(F)*Fmin/max(F))+1;
  plot(real(Xfer(imin:imax)),imag(Xfer(imin:imax)))
  xlabel('Real')
  ylabel('Imaginary')
  grid on
  zoom on
 elseif FLAG==5
  subplot(1,1,1)
  mag=20*log10(abs(Xfer));
  plot(F,mag)
  xlabel('Frequency (Hz)')
  ylabel('Mag (dB)')
  grid on
  axis([Fmin Fmax min(mag) max(mag)])
  zoom on
 elseif FLAG==6
  subplot(1,1,1)
  phase=unwrap(angle(Xfer))*180/pi;
  plot(F,phase)
  xlabel('Frequency (Hz)')
  ylabel('Phase (deg)')
  grid on
  phmin_max=[floor(min(phase)/45)*45 ceil(max(phase)/45)*45];
  axis([Fmin Fmax  phmin_max(1) phmin_max(2)])
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
  zoom on
 elseif FLAG==7
  subplot(1,1,1)
  plot(F,real(Xfer))
  xlabel('Frequency (Hz)')
  ylabel('Real')
  grid on
  axis([Fmin Fmax min(real(Xfer)) max(real(Xfer))])
  zoom on
 elseif FLAG==8
  subplot(1,1,1)
  plot(F,imag(Xfer))
  xlabel('Frequency (Hz)')
  ylabel('Imaginary')
  grid on
  axis([Fmin Fmax min(imag(Xfer)) max(imag(Xfer))])
  zoom on
 else
  subplot(2,1,1)
  mag=20*log10(abs(Xfer));
  plot(F,mag)
  xlabel('Frequency (Hz)')
  ylabel('Mag (dB)')
  grid on
  axis([Fmin Fmax min(mag) max(mag)])
  zoom on
  subplot(2,1,2)
  phase=unwrap(angle(Xfer))*180/pi;
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
