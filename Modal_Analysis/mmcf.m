function [z,nf,u]=mmcf(f,H,Fmin,Fmax,nmodes)
%[z,nf,u]=mmcf(f,TF,Fmin,Fmax,nmodes) Curve fit to MDOF FRF using polyreference.
% f is the frequency vector in Hz. 
%
% H is the complex frequency response functions, each FRF in a column.
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
% Hacked 12/18/01 to no avail. Still doesn't work.



inlow=1;
inhigh=length(f);
if nargin==4 |nargin==5
inlow=floor(length(f)*(Fmin-min(f))/(max(f)-min(f)))+1;
inhigh=ceil(length(f)*(Fmax-min(f))/(max(f)-min(f)))+1;
elseif nargin==3
	nmodes=Fmin;
end

if f(inlow)==0
 inlow=1;
end
size(f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we make No x Ni blocks at the same frequencies when we have 
% MIMO data. 
Ni=size(H,2);
No=size(H,3);
if size(H,3)~=1
	f=f(inlow:inhigh,:);
	H=H(inlow:inhigh,:,:);

	Hold=H;
	fold=f;
	M=size(H,1);
% If No>Ni, we need to fix this. This assumes reciprocity holds.
	if No>Ni
		H=zeros(Hold);
		for l=1:Ni
			H(:,:,l)=squeeze(H(:,l,:));
		end
	end
	H=zeros(M*No,Ni);
	for l=1:No
		H(l:No:No*M,:)=Hold(:,:,l);
		f(l:No:No*M,1)=fold(:,1);
	end
else
	f=f(inlow:inhigh,:);
	
	H=H(inlow:inhigh,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tfplot(f(3:3:length(f),:),H(3:3:length(f),:))
% break
% Reduce to nmodes FRF for finding poles:
max(f),min(f)
tfplot(f(3:3:2400),H(3:3:2400,1))
size(H),size(f)
return
R=(H');
[U,S,V]=svd(R);

svals=diag(S)
if nargin==2|nargin==4
	lenS=length(svals);
	countmodes=svals(2:lenS)>svals(1:lenS-1)*.05;
	nmodes=sum(countmodes)+1;
	disp([num2str(nmodes) ' modes found.'])
end
numfreq=nmodes;
T=U(:,1:nmodes);
Hp=(T')*R;
R=(Hp');
H;
%disp('R=H')
%R=H;
%tfplot(f,R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From here on, the reduced FRFs are stored in R. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we make the negative frequency parts of the FRF (and frequency vector)
% This allows us to obtain both poles in each pair.
M=length(f)/No; % This change when we zoomed in frequency
R=[conj(R(M*No:-1:1,:,:));R];
w=2*pi*[-f(M*No:-1:1);f];
M=length(w)/No;
%tfplot(abs(w)/2/pi,R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to obtain the companion matrix, we need to form the 
% H omega matrix. We then solve for the alpha-beta matrix, then find 
% eigenvalues and eigenvectors

% order of the method being used
% 1 is Craig (Multi Reference Frequency Domain)
% % Polyreference Frequency Domain
m=2; 
diagomega=spdiags(w*sqrt(-1),0,M*No,M*No);
INo=spdiags(ones(No,M*No),[-No*(M-1):No:0],No*M,No);

HO=zeros(M*No,Ni*m+No*(m+2));
%HO=sparse(HO);
% Since we need one more beta column, and its power is zero, we 
% put it in up front.

% We should fill in the one lowest H, and the 3 lowest Omegas,
% since they are always there (in the generic first order case)
% and loop only for the other terms.
%
% The higher Omega parts (beta) are not shown in Allemang's
% paper, but they represent residual terms.
%
%tfplot(f,R)
R;
HO(:,1:Ni)=R;
%pause
HO(:,m*Ni+(1:No))=-INo;
HO(:,m*Ni+No+(1:No))=-diagomega*INo;
HO(:,m*Ni+2*No+(1:No))=-diagomega.^2*INo;


for l=2:m
	HO(:,l*Ni+(1-Ni:0))=diagomega.^(l-1)*R;
	HO(:,m*Ni+(l+1)*No+(1:No))=-diagomega.^(l+1)*INo;
end
% We still need the right hand side of the equation
RHS=-diagomega.^m*R;
%disp('HO')
%HO
%RHS
% Now the alpha-beta matrix comes from the least squares solution
alphabeta=HO\RHS;
%size(alphabeta)
%(HO*alphabeta-RHS)./RHS
%break
% Here we need to assemble the state space (companion) matrix.
if m==1
	C=-alphabeta(1:Ni,1:Ni);	
else
	C=spdiags(ones(Ni*(m-1),1),-Ni,Ni*m,Ni*m);
%	full(C)
	
    for im=1:m
		C(1:Ni,(m-im)*Ni+(1:Ni));
		alphabeta((1:Ni)+Ni*(im-1),:);
		C(1:Ni,(m-im)*Ni+(1:Ni))=-alphabeta((1:Ni)+Ni*(im-1),:);
	end
end
%full(C)
%%%%%%%%%%%%%(^%(&*^(^%^%*^%*^%()))*%^^)((((((((((((((((((((((((*^%&(^%*&^%&*^%&*%^*%^$

% Here we solve for eigenvalues and eigenvectors

damp(full(C))
% Here we find the modal participation factors

disp('What happened here?')

%break
% Here we translate the mode shapes and modal participation factors back into real space.

%               junk below here
n=1;
N=2*n;
N=2;
%[x,y]=meshgrid(0:N,R);
%[x,w2d]=meshgrid(0:N,w);
%disp('make c')
tic
%cold=diag(-w.^(N))*R;
%size(c)
%toc
%disp('alternative make c')
%numfreq=
tic
numfreq=1023;
diagomega=spdiags(w,0,numfreq*2,numfreq*2);
size(diagomega)
size(R)
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
%%%%%%%%%%%%%%    junk above here
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

H3=[H2;H];

%plot(unwrap(angle(H3)))
a=XoF\H3;


%size(w)
%plot(abs(w)),pause
%tfplot(f,H)
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
	%plot([abs(XoF) abs(H)])
	%2*ll
	%plot(abs(w(ll:2*ll)))
	%pause
	%break

	Fmin=min(f);
	Fmax=max(f);
	phase=unwrap(angle(H(:,mm)),pi*1.5)*180/pi;
	phase2=unwrap(angle(XoF),pi*1.5)*180/pi;size(phase);
	%size(XoF)
	subplot(2,1,1)
	plot(f,20*log10(abs(XoF)),f,20*log10(abs(H(:,mm))))

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
