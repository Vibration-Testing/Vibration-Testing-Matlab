function [z,nf,a]=sdofcf(f,TF)
%[z,nf,a]=sdof(f,TF) Curve fit to SDOF FRF.
% f is the frequency vector in Hz. It does not have to 
%    start at 0 Hz.
% TF is the complex transfer function.
% z and nf are the damping ratio and natural frequency (Hz)
% a is the product of the residues of the coordinates the 
% transfer function is between. (For example, in the example 
% below, a1 times a2 is returned. See equation 7.42)
% Only one peak may exist in the segment of the FRF passed to 
% sdofcf. No zeros may exist withing this segment. If so, 
% curve fitting becomes unreliable. 
%
% EXAMPLE:
% M=eye(2);
% K=[2 -1;-1 2];
% C=.01*K;
% [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,2,linspace(0,.5,1024));
% figure(1)
% n=250;
% f2=Freq((1:n)+450);
% R2=Recep((1:n)+450);
% R2=R2+.1*randn(n,1)+.1*randn(n,1)*i;% Poorly Simulated Noise
% [z,nf,a]=sdofcf(f2,R2)
%
% Note that by changing the parts of Freq and Recep used
% We can curve fit to other modes.

% Copyright Joseph C. Slater, 10/8/99
% Updated 11/8/99 to improve robustness
% Updated 1/1/00 to improve robustness
% Updated 4/1/01 from vtb7_4 to linear curve fitting.
%tfplot(f,TF)
%pause
R=TF;
[y,in]=max(abs(TF));

ll=length(f);
w=f*2*pi*sqrt(-1);
w2=w*0;
R3=R*0;
for i=1:ll
	R3(i)=conj(R(ll+1-i));
	w2(i)=conj(w(ll+1-i));
end
w=[w2;w];
R=[R3;R];

c=-w.^2.*R;
a0=w.^0.*R;
a1=w.^1.*R;
a2=R*0-1;
a=[a0 a1 a2];
%a(1:10,:)
b=a\c;
rs=roots([1 b(2) b(1)]);

omega=abs(rs(1));
z=-real(rs(1))/abs(rs(1));
nf=omega/2/pi;

XoF=1./(w.^2+2*z*w*omega+omega^2);

a=1;

a=XoF\R;
%tfplot(f,TF)
%tfplot(f,[a*Rest,R])
XoF=a*XoF(ll+1:2*ll);


%break
Fmin=min(f);
	  Fmax=max(f);
	  phase=unwrap(angle(TF))*180/pi;
	  phase2=unwrap(angle(XoF))*180/pi;size(phase);
		%size(XoF)
	  subplot(2,1,1)
	  plot(f,20*log10(abs(XoF)),f,20*log10(abs(TF)))
	  as=axis;
	  zoom on
	  legend('Identified FRF','Experimental FRF',0)
	  min(f);

		axis([Fmin Fmax as(3) as(4)])
	  xlabel('Frequency (Hz)')
	  ylabel('Mag (dB)')
	  grid on
	%  Fmin,Fmax,min(mag),max(mag)
	%  axis([Fmin Fmax minmag maxmag])
	  while phase2(in)>50
		  phase2=phase2-360;
	  end
	  phased=phase2(in)-phase(in);
	  phase=phase+round(phased/360)*360;
	  phmin_max=[floor(min(min([phase;phase2]))/45)*45 ceil(max(max([phase;phase2]))/45)*45];
	  subplot(2,1,2)
	  plot(f,phase2,f,phase)
	  xlabel('Frequency (Hz)')
	  ylabel('Phase (deg)')
	  legend('Identified FRF','Experimental FRF',0)
	
	  grid on
	  axis([Fmin Fmax  phmin_max(1) phmin_max(2)])
	  gridmin_max=round(phmin_max/90)*90;
	  set(gca,'YTick',gridmin_max(1):22.5:gridmin_max(2))
	  zoom on
	
