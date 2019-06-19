function [z,nf,u]=mmcf(f,TF,Fmin,Fmax,nmodes)
%[z,nf,u]=mmcf(f,TF,Fmin,Fmax,nmodes) Curve fit to MDOF FRF.
% f is the frequency vector in Hz. 
%
% TF is the complex transfer functions, each FRF in a column.
% z and nf are the damping ratio and natural frequency (Hz)
% u is the mode shape.
%
% If Fmin and Fmax are not entered, the min and max values in
% f are used.
%
% If nmodes is not included, it is estimated automatically.
%
% EXAMPLE:
% M=eye(2);
% K=[2 -1;-1 2];
% C=.01*K;
% [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,2,linspace(0,.5,1024));
% R1=Recep+.1*randn(1024,1)+.1*randn(1024,1)*i;% Poorly Simulated Noise
% [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,1,linspace(0,.5,1024));
% R2=Recep+.1*randn(1024,1)+.1*randn(1024,1)*i;% Poorly Simulated Noise
% [z,nf,u]=mmcf(Freq,[R1 R2])
%
% Note that by changing the parts of Freq and Recep used
% We can curve fit to other modes.

% Copyright Joseph C. Slater, 10/8/99
% Updated 11/8/99 to improve robustness
% Updated 1/1/00 to improve robustness
% Updated 4/1/01 from vtb7_4 to linear curve fitting.
% Updated 5/1/01 to start adding MDOF.
% Updated 5/15/01 to include frequency range for curve fitting
% Modified from MDOFCF 7/1/01 to include multiple modes

disp('This is beta software written 02-Jun-2001')
disp('This code may not be used by anyone not authorized by J. Slater')
disp('This code will expire 01-Jan-2002')

if datenum('10-Dec-2001')<datenum(date)
	warndlg('contact joseph.slater@wright.edu for an update.','mdofcf.p will expire 01-Jan-2002')	
end


if datenum('01-Jan-2002')<datenum(date)
	warndlg('contact joseph.slater@wright.edu for an update.','mdofcf.p has expired')	
end

% Cut above here for p-code
if datenum('01-Jan-2002')<datenum(date)
	delete mdofcd.p
end

inlow=1;
inhigh=length(f);
if nargin==4 |nargin==5
inlow=floor(length(f)*(Fmin-min(f))/(max(f)-min(f)))+1;
inhigh=ceil(length(f)*(Fmax-min(f))/(max(f)-min(f)))+1;
elseif nargin==3
	nmodes=Fmin;
end

if f(inlow)==0
 inlow=2;
end
if size(TF,3)~=1
	f=f(inlow:inhigh,:);
	TF=TF(inlow:inhigh,:,:);

	TFold=TF;
	fold=f;
	s1=size(TF,1);
	s2=size(TF,2);
	s3=size(TF,3);
	TF=zeros(s1*s3,s2);
	for l=1:s3
		TF(l:s3:s3*s1,:)=TFold(:,:,l);
		f(l:s3:s3*s1,1)=fold(:,1);
    end
else
	f=f(inlow:inhigh,:);
	TF=TF(inlow:inhigh,:);
end
%tfplot(f,TF)

% Reduce to nmodes FRF for finding poles:
R=(TF');
[U,S,V]=svd(R);


svals=diag(S);
if nargin==2|nargin==4
	lenS=length(svals);
	countmodes=svals(2:lenS)>svals(1:lenS-1)*.05;
	nmodes=sum(countmodes)+1;
	disp([num2str(nmodes) ' modes found.'])
end
T=U(:,1:nmodes);
Hp=(T')*R;
R=(Hp');
tfplot(f,R)
nmodes
break
%tfplot(f,R)
%disp('This line needs to be removed or corrected')
%R=(TF);

% End reduction (whole thing is still stored in TF)
w=f*2*pi*sqrt(-1);
numfreq=length(w);
%This replaces the next section, which is rather slow.
%tic
w3=zeros(size(w,1)*2,1);
R4=zeros(size(w,1)*2,size(R,2));
TF3=zeros(size(w,1)*2,size(TF,2));
w3=[conj(w(numfreq:-1:1,:));w];
R4=[conj(R(numfreq:-1:1,:));R];
TF3=[conj(TF(numfreq:-1:1,:));TF];
%toc

% Now we make the negative frequency parts of the FRF (and frequency vector)
ll=length(f);
tic
w3=zeros(size(w,1)*2,1);
w2=w*0;
R3=R*0;

for i=1:ll
	R3(i,:)=conj(R(ll+1-i,:));
	w2(i,:)=conj(w(ll+1-i,:));
	TF2(i,:)=conj(TF(ll+1-i,:));
end
w=[w2;w];
R=[R3;R];
%toc
%tfplot(imag(w)/2/pi,[R R4]);

n=1;
N=2*n;
N=2;
%[x,y]=meshgrid(0:N,R);
%[x,w2d]=meshgrid(0:N,w);
%disp('make c')
%tic
%cold=diag(-w.^(N))*R;
%size(c)
%toc
%disp('alternative make c')
%tic
diagomega=spdiags(w,0,numfreq*2,numfreq*2);
c=-diagomega.^2*R;
%size(c)
%[w full(diag(diagomega))]
%plot(real(diag(diagomega)),real(w)),pause,break


%toc
% The following line was from mdofcf
%aa=[w2d(:,1+(0:N-1)).^x(:,1+(0:N-1)).*y(:,1+(0:N-1)) -w2d(:,1+(0:N)).^x(:,1+(0:N))];

%Replacement line for previous line to handle a matrix, the columns of which are
%individual FRFs. Tested first half
% Right half unrewritten (from -w2d on).

%The following line was an intermediate step.
%aa2=[ R       diag(w.^1)*R      -w2d(:,1+(0:N)).^x(:,1+(0:N))];

% Rewrite of above line.
sizeR=size(R,2);
%[junk,wspread]=meshgrid(size(R,2),w);
tic
%disp('make aa')
aa=[ R       diagomega*R    -w.^0  -w.^1 -w.^2];
%toc
%more on
%aa2-aa3
%aa2=aa3;toc
%aa=aa2;toc
%disp('get b'),tic
b=aa\c;%toc

%[1 b((N-1:-1:0)+1)'];


%rs=roots([1 b((N-1:-1:0)+1)'])
tic
%disp('solve eigenvalue problem')
b=conj(b');
[xx,ee]=polyeig(b(:,1:sizeR),b(:,sizeR+1:sizeR*2),eye(sizeR));
rs=ee;
%toc
% We want to sort eigenvalues in increasing order, prioritizing
% the imaginary parts. Matlab sorts by magnitude.
% Here we take only every other root since the other is the conjugate
% and therefore redundant.
rs(1:2:length(rs));
[yrs,irs]=sort(abs(imag(rs)));

rs=rs(irs(1:2:sizeR*2));

omega=abs(rs);
z=-real(rs)./abs(rs);
nf=omega/2/pi;
%toc


% This section now gets repeated for each FRF. Still needs to be done for MDOF, to make a matrix a.
%rs)
%XoF=[1./(w-rs(1)) 1./(w-rs(2)) 1./(w.^0) 1./w.^2];%;1./w2d(:,[1 3])];

%disp('OK, just expand this to multiple poles')
[rsmesh,wmesh]=meshgrid([rs; conj(rs)],w);
XoF1=1./(wmesh-rsmesh);

%XoF1=[1./(w-rs(1)) 1./(w-rs(2)) 1./(w-conj(rs(1))) 1./(w-conj(rs(2)))];
%XoF2=1./(w+rs(1));
XoF2=1./(w.^0);
XoF3=1./w.^2;
XoF=[XoF1 XoF2 XoF3];

%for i=1:

TF3=[TF2;TF];

%plot(unwrap(angle(TF3)))
a=XoF\TF3;


%size(w)
%plot(abs(w)),pause
%tfplot(f,TF)
%tfplot(f,[a*Rest,R])
%size(XoF)
u=a(1:sizeR,:)';
%u=u/u(1);
%u=u/norm(u);
u=u/diag(exp(angle(u(1,:))*sqrt(-1)));
R=XoF;
for mm=1:size(a,2)

	if mm>1
		pause
	end
clf
	figure(gcf)
	XoF=R((ll+1:2*ll),:)*a(:,mm);
	%plot([abs(XoF) abs(TF)])
	%2*ll
	%plot(abs(w(ll:2*ll)))
	%pause
	%break

	Fmin=min(f);
	Fmax=max(f);
	phase=unwrap(angle(TF(:,mm)),pi*1.5)*180/pi;
	phase2=unwrap(angle(XoF),pi*1.5)*180/pi;size(phase);
	%size(XoF)
	subplot(2,1,1)
	plot(f,20*log10(abs(XoF)),f,20*log10(abs(TF(:,mm))))

	as=axis;
	legend('Identified FRF','Experimental FRF',0)
	min(f);
	axis([Fmin Fmax as(3) as(4)])
	title(['Frequency Response Function ' num2str(mm) ' Fit'])
	xlabel('Frequency (Hz)')
	ylabel('Mag (dB)')

	grid on
	zoom on

	drawnow

	%  Fmin,Fmax,min(mag),max(mag)
	%  axis([Fmin Fmax minmag maxmag])
	%pause
	in=floor(length(f)/2);
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
	
	axis([Fmin Fmax  phmin_max(1) phmin_max(2)])
	gridmin_max=round(phmin_max/90)*90;
	set(gca,'YTick',gridmin_max(1):22.5:gridmin_max(2))
	grid on
	zoom on
	drawnow
	figure(gcf)
end
