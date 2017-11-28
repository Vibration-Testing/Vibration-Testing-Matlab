function [freqout,tfout]=tfest(x,f,dt,n,options)
%TFEST Estimates Transfer Function
% [Freq,Txf] = TFEST(X,F,DT,N,OPTIONS) estimates the 
% Frequency Response Function (FRF) between X and F 
% (X/F) using either a Fast Fourier Transform when 
% possible, or a Discrete Fourier Transform otherwise.
% TFEST calculates the H1 FRF (Sxf/Sff) by default.
% Freq is the frequency vector in Hertz and Txf is the
% Transfer Function in complex form. N is the number 
% of points used in the Fourier transform. The default 
% value is the length of the vectors if it is a power 
% of 2 or the next power of 2 after the length of the 
% vectors. DT is the time step of the sampled data. 
% If DT is replaced by the time vector, DT will be extracted 
% from the time vector using DT = T(2) - T(1). If X and 
% F are matrices, TFEST  will find the frequency response
% function for each column and average the results unless 
% OPTIONS{1} is set to 'no'. 
% The value of OPTIONS{2} defines the type of FRF estimation.
% OPTIONS{2} = 1 is the default, resulting in use of the H1 
% FRF algorithm (Sxf/Sff) minimizing noise on the output. If
% OPTIONS{2} = 2, the H2 FRF algorithm (Sxx/Sff) is used,
% minimizing noise on the input. If OPTIONS{2} = 3, the Hv 
% algorithm is used, minimizing noise on both the input and 
% the output. 
% 
% NOTE: Using the Hv algorithm will take significantly longer
% to run.
% NOTE: OPTIONS is a cell array. Please be sure to use curly brackets "{}"
% as shown if you choose to use this option.
% N and OPTIONS are optional. Either can be
% be left out.
%
% TFEST(X,Y,DT,N,OPTIONS) plots the Frequency Response Function 
% if there are no output arguments.  Click in the region of 
% interest to zoom in.  Each click will double the size of the 
% plot.  Double click to return to full scale.
%
% See also COH, ASD, CRSD, and TFPLOT.

% Copyright (c) 1994 by Joseph C. Slater
% Modifications:
% --------------
% 7/6/00: Changed default FRF calculation from H2 to H1
%         Added H1, H2, and Hv options.
sy=size(f);
if nargin==3
  n=2^nextpow2(sy(1));
  ave='yes';
  frftype=1; %H1
elseif nargin==4
  if strcmp(n(1),'no')|strcmp(n(1),'yes')
    ave=n{1};
    if length(n)==1
      frftype=1;
    else
      frftype=n{2};
    end
    n=2^nextpow2(sy(1));
  else
    ave='yes';
    frftype=1;
  end
else
  frftype=options{2};
  ave=options{1};
end

if isempty(n)
  n=2^nextpow2(sy(1));
end


if length(dt)~=1
  dt=dt(2)-dt(1);
end


if frftype==1
  [freq,Pff]=asd(f,dt,n);
  [freq,Pxf]=crsd(x,f,dt,n);
  tfunc=Pxf./Pff;
  disp('H1')
elseif frftype==2
  [freq,Pxx]=asd(x,dt,n);
  [freq,Pxf]=crsd(x,f,dt,n);
  tfunc=Pxx./conj(Pxf);
  disp('H2')
elseif frftype==3
  [freq,Pxx]=asd(x,dt,n);
  [freq,Pff]=asd(f,dt,n);
  [freq,Pxf]=crsd(x,f,dt,n);
  disp('Hv')
  for i=1:size(Pxx,1)
    for j=1:size(Pxx,2)
      frfm=[Pff(i,j) conj(Pxf(i,j));Pxf(i,j) Pxx(i,j)]
      [v,d]=eig(frfm)
      [y,yi]=sort(diag(d));
      tfunc(i,j)=-v(1,yi(1))/v(2, yi(1));%,pause		
    end
  end
end

tfunc=conj(tfunc);

sfft=size(Pxf);
if sfft(2)~=1 & ~strcmp(ave,'no')
  tfunc=mean(tfunc')';
end

% If no left hand arguments then plot results
if nargout==0
  subplot(211)
  plot(freq,20*log10(abs(tfunc)))
  
  xlabel('Frequency (Hz)')
  ylabel('Mag (dB)')
  grid
  zoom on
  subplot(212)
  phase=unwrap(angle(tfunc))*180/pi;
  %phase=angle(tfunc)*180/pi;
  plot(freq,phase)
  
  xlabel('Frequency (Hz)')
  ylabel('Phase (deg)')
  grid
  phmin_max=[floor(min(phase)/45)*45 ceil(max(phase)/45)*45];
  set(gca,'YLim',phmin_max)
  gridmin_max=round(phmin_max/90)*90;
  set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))

  set(gca,'GridLineStyle',':')
  set(gca,'YTickLabel',gridmin_max(1):90:gridmin_max(2))
  zoom on
  return
end

freqout=freq;
tfout=tfunc;
