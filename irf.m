function [tout,impout]=irf(frf,dF,n,ave)
%IRF Impulse Response Function estimate of the Frequency 
% response function FRF.
% [T,I] = IRF(FRF,DF,N,AVE) estimates the Impulse Response 
% Function by taking the inverse DFT of the Frequency 
% Response Function. N is the number of points to be used 
% in the Inverse Fourier Transform. The default is N=
% length(FRF). If DF is the frequency vector in Hz, the 
% frequency spacing between FRF points is extracted as 
% F(2) - F(1). If FRF is a matrix, IRF will find the 
% Impulse Response Function for each column and average 
% the results unless AVE is set to 'noave'. N and AVE 
% are optional. Either one may be left out.
%
% IRF(FRF,DF,N,AVE) plots the Impulse Response Function if 
% there are no output arguments. 
%
% See Also TFEST

% Copyright (c) 1994 by Joseph C. Slater

sfrf=size(frf);

if nargin==2
  ave='yes';
  n=sfrf(1);
 elseif nargin==3
  if strcmp(n,'noave')
    ave=n;
    n=length(frf);
   else
    ave='yes';
  end
end
if n==[]
  n=2^nextpow2(sfrf(1));
end
 
if length(dF)>1
 dF=dF(2)-dF(1);
end
lfrf=sfrf(1);

frfmir=[frf ; zeros(1,sfrf(2)); conj(frf(length(frf):-1:2,:))];

impres=ifft(frfmir,n*2);

lenimpres=length(impres);
dt=1/dF;
t=(0:lenimpres-1)'/(lenimpres-2)*dt;


if sfrf(2)~=1 & ~strcmp(ave,'noave')
 impres=mean(impres')';
end

if nargout==0
 plot(t,real(impres))

 title('Impulse Response')
 xlabel('Time')
 ylabel('Impulse Response')
 grid
 zoom on
 if norm(imag(impres))>norm(real(impres))/100
   warndlg(['Frequency response data is faulty. Impulse  ';...
            'response is complex. Please check your input';...
    		   'data. Continuing with calculations.         ']);
 end
% impres;
 return
end

if norm(imag(impres))>norm(real(impres))/100
  warndlg(['Frequency response data is faulty. Impulse  ';...
           'response is complex. Please check your input';...
  		    'data. Continuing with calculations.         ']);
end
impres;
impres=real(impres);

tout=t;
impout=impres;
