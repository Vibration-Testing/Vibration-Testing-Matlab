function [z,nf,amp]=vtb7_4_1(f,TF1,TF2,TF3,TF4,TF5,TF6,b)
%
%[z,nf,amp]=vtb7_4_1(f,TF1,TF2,TF3,TF4,TF5,TF6) Curve fit to MDOF FRF.
% f is the frequency vector in Hz. It does not have to 
% start at 0 Hz.
% TF1 is the complex transfer function for node 1.
% TF2 is the complex transfer function for node 2.
% TF3 is the complex transfer function for node 3.
% TF4 is the complex transfer function for node 4.
% TF5 is the complex transfer function for node 5.
% TF6 is the complex transfer function for node 6.
% z and nf are the damping ratio and natural frequency (Hz)
% amp is the amplitude vector giving the amplitude at each 
% node corresponding to the natural frequency.
%
% STOP! To be able to use this code, you need to have the complex
% transfer function data corresponding to all the six nodes. For this
% purpose run the code 'loaddata'.

%##################################################################### 
global XoF
global sz
if nargin==7
    
    
  
    TF=[TF1,TF2,TF3,TF4,TF5,TF6];
    sz=max(size(TF(1,:)));
    lf=length(f);
    z=0.0005;	
    xx=zeros(sz+2+3*sz,1); 
     
    for J=1:sz		
    TFN=TF(:,J);
    [y,in]=max(abs(TFN));
    f(in);
    a0=abs(TFN(1))*(2*pi*f(in))^2;
    a0=-sign(imag(TFN(in)))*abs(TFN(in))*2*z*(2*pi*f(in))^2;
    xx(J,1)=a0;
    end
    
    xx(sz+1,1)=z;
    xx(sz+2,1)=2*pi*f(in);
    x=xx;
 
   
	
    if in-3<1|in+2>length(f)
		disp('The peak response must be near the middle of the data')
		disp('Please center your peak and try again')
		break
    end
	


x=fmins('vtb7_4_1',x,[],[],f(in-3:in+2),TF1(in-3:in+2),TF2(in-3:in+2),TF3(in-3:in+2),TF4(in-3:in+2),TF5(in-3:in+2),TF6(in-3:in+2));
	x;
	
	cost=vtb7_4_1(x,f,TF1,TF2,TF3,TF4,TF5,TF6);
	x=fmins('vtb7_4_1',x,[],[],f,TF1,TF2,TF3,TF4,TF5,TF6);
	x;
	
	cost=vtb7_4_1(x,f,TF1,TF2,TF3,TF4,TF5,TF6);
	x=fmins('vtb7_4_1',x,[],[],f,TF1,TF2,TF3,TF4,TF5,TF6);
	x;
	
	cost=vtb7_4_1(x,f,TF1,TF2,TF3,TF4,TF5,TF6);
	x=fmins('vtb7_4_1',x,[],[],f,TF1,TF2,TF3,TF4,TF5,TF6);
	x;
	
	cost=vtb7_4_1(x,f,TF1,TF2,TF3,TF4,TF5,TF6);
	x=fmins('vtb7_4_1',x,[],[],f,TF1,TF2,TF3,TF4,TF5,TF6);
	x; 
	
	cost=vtb7_4_1(x,f,TF1,TF2,TF3,TF4,TF5,TF6);
	x=fmins('vtb7_4_1',x,[],[],f,TF1,TF2,TF3,TF4,TF5,TF6);
	x;
	
	z=x(sz+1);
	om=x(sz+2);
	nf=om/2/pi;	
	amp=x(1:sz); 
	plot(amp);
        grid
	title('Mode Shape Plot')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
else
	x=f;
	f=TF1;
	TF1=TF2;
	TF2=TF3;
 	TF3=TF4;
	TF4=TF5;
	TF5=TF6;
	TF6=b;
	
	w2=f*2*pi;
	lx=length(x);
	x(sz+1)=abs(x(sz+1));
	x(sz+2)=abs(x(sz+2));
	
	
	xfer=[TF1,TF2,TF3,TF4,TF5,TF6];
	t=1;
	error1=0;
	for K=sz+3:3:lx
	XoF=x(K)+x(K+1)*i*w2-x(K+2)*w2.^2;
	XoF=XoF+x(t)./(-w2.^2+2*x(sz+1)*w2*i*x(sz+2)+x(sz+2)^2);
	error=norm(XoF-xfer(:,t));
	error1=error1+error;
	t=t+1;
	end
	sla74=error1;
	z=sla74;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
