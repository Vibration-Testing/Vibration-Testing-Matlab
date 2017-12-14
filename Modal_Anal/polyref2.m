function [z,nf,poles,u]=polyref2(f,H,Fmin,Fmax,nmodes,options)
%[z,nf,poles,u]=mmcf(f,H,Fmin,Fmax,nmodes) Curve fit to MDOF FRF using polyreference.
% f is the frequency vector in Hz. 
% H is the matrix of FRFs, in blocks of p outputs by o inputs, stacked
%   vertically by increasing frequency for each block.
%
% If Fmin and Fmax are not entered, the min and max values in
% f are used.
% nmodes is the number of modes to presume
% If nmodes is left out, svd of the data is performed, and cool stuff
% happens. Eventually
% z and nf are the damping ratio and natural frequency (Hz)
% p are the poles.
% u are the mass normalized mode shapes,
% options=[viewsvd reviewfit debug]
%    viewsvd=1 plots svds and prompts the user for a number of modes for
%    viewsvd=2 attempts to chose number of modes automatically
%    reviewfit=1 plots the reduced FRFs against the curve fits, 1 at a time
%    debug=1 causes a whole bunch of diagnostics to be displayed
%
% EXAMPLE:
%
%  E=7.31e10;
%  I=1/12*.03*.015^3;
%  rho=2747;
%  A=.015*.03;
%  L=0.5;
%  [f,h1]=vtb6_5(.1,.1,0,400,[E I rho A L],2,0.005,200);
%  [f,h2]=vtb6_5(.1,.15,0,400,[E I rho A L],2,0.005,200);
%  [f,h3]=vtb6_5(.1,.2,0,400,[E I rho A L],2,0.005,200);
%  [f,h4]=vtb6_5(.1,.25,0,400,[E I rho A L],2,0.005,200);
%  H=[h1 h2 h3 h4];
%  tfplot(f,H)
%  [z,nf,poles,u]=polyref2([f ],H,2)
%
% EXAMPLE:
% M=eye(2);
% K=[2 -1;-1 2];
% C=.01*K;
% [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,2,linspace(0,.5,1024));
% R1=Recep;%+.1*randn(1024,1)+.1*randn(1024,1)*i;% Poorly Simulated Noise
% [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,1,linspace(0,.5,1024));
% R2=Recep;%+.1*randn(1024,1)+.1*randn(1024,1)*i;% Poorly Simulated Noise
% [z,nf,poles,u]=polyref2(Freq,[R1 R2],2)

% Copyright Joseph C. Slater, 2/24/12

global freqdebug
freqdebug=0;
reviewfit=0;
viewsvd=0;
if nargin==6
    if length(options)>1
        reviewfit=options(2);
    end
    if length(options)>2
        freqdebug=options(3);
    end
    viewsvd=options(1);
end

z=[];
nf=[];
u=[];
p=size(H,1)/length(f);
o=size(H,2);
if length(f)~=size(H,1)
    if p~=floor(p)
        disp(['You cannot have ' num2str(size(H,2)/length(f)) ' outputs'])
        return
    end    
end

if nargin==2
    nmodes =size(H,2);
elseif nargin==3
    nmodes=Fmin;
elseif nargin==4
    nmodes =size(H,2);
    inlow=max([floor(length(f)*(Fmin-min(f))/(max(f)-min(f)))+1 1]);
    inhigh=min([ceil(length(f)*(Fmax-min(f))/(max(f)-min(f))) length(f)]);
    f=f(inlow:inhigh);
    H=H((1+(inlow-1)*p):(inhigh*p),:);
elseif nargin==5|nargin==6
    Fmin=max([Fmin min(f)]);
    Fmax=min([Fmax max(f)]);
    inlow=max([floor(length(f)*(Fmin-min(f))/(max(f)-min(f)))+1 1])
    inhigh=min([ceil(length(f)*(Fmax-min(f))/(max(f)-min(f))) length(f)])
    size(f)
    
    f=f(inlow:inhigh);
    H=H((1+(inlow-1)*p):(inhigh*p),:);
end

if freqdebug==1, disp('Appending FRF'),end

f=[-f ;f];
H=[conj(H);H];


if p>o
    if freqdebug==1, disp('Re-sorting H'),end
    
    H=Vscript(H,f*2*pi,0);
    disp('p>o')
end

Horiginal=H;
% Here we reduce the FRFs to the number of modes the user thinks we have.
if freqdebug==1, disp('Reducing FRFs'),end

[U,S,~]=svd(H,'econ');
lf=length(f);
percentcontent=0;
ii=0;
ds=diag(S);
while percentcontent<.75&&ii<length(ds)
    ii=ii+1;
    percentcontent=sum(ds(1:ii))/sum(ds);
end

    
while viewsvd==1
    
    clf
    subplot(2,1,1)
    plot(diag(S),'*')
    axis([0 length(ds)+1 0 max(ds)*1.1])
    grid on
    title('Singular values (strength of individual components)')
    subplot(2,1,2)
    plot(f(lf/2+1:lf),20*log10(abs(U(lf/2+1:lf,1))))
    grid on
    xlabel('Frequency (Hz)')
    ylabel('FRF (dB)')
    title('FRF of Principle FRF')
    figure(gcf)
    nmodes2=input(['How many modes do you think there are? (0 to abort, return for auto ' num2str(ii) '.) ']);
    if nmodes2==0
        disp('Aborting curve fit')
        z=[];nf=[];poles=[];u=[];
        return
    elseif nmodes2==-1
        input('Enter frequency range as [Fmin Fmax]')
        viewsvd=1;
    elseif isempty(nmodes2)
        nmodes=ii;
    else 
        nmodes=nmodes2;
    end
    viewsvd=0;
end

if viewsvd==2
    nmodes=ii;
end



H=U(:,1:nmodes);


%for iii=1:nmodes
%    tfplot(f,H(:,iii))

%    pause
%end



if freqdebug==1, disp('Plot and look at reduced FRFs'),end

if length(f)==size(H,1)
    if freqdebug==1, disp('Re-order H'),end
    H=Vscript(H,f*2*pi,0);
end

if freqdebug==1, disp('Make Vscript0_1'),end
Vscript0_1=Vscript(H,f*2*pi,[0 1]);

if freqdebug==1, disp('Make Uscript0_1'),end
Uscript0_1=Uscript(size(H,2),f*2*pi,[0 1 2 ]);

if freqdebug==1, disp('Make Vscript2'),end
Vscript2=Vscript(H,f*2*pi,2);

if freqdebug==1, disp('Performing Gauss elimination'),end
AA=[Vscript0_1 -Uscript0_1]\(-Vscript2);
if freqdebug==1, disp('Finished Gauss elimination'),end
%(-Vscript(H,f*2*pi,2))-[Vscript(H,f*2*pi,[0 1]) -Uscript(size(H,2),f*2*pi,[0 1 2 ])]*AA
sA=size(AA,2);
sA1=size(AA,1);
sB=(sA1-2*sA)/3;
AA;
A0=AA(1:sA,1:sA).';
A1=AA(sA+(1:sA),1:sA).';
B0=-AA(sA*2+(1:sB),1:sA).';
B1=-AA(sA*2+sB+(1:sB),1:sA).';
B2=-AA(sA*2+sB*2+(1:sB),1:sA).';

%sqrt(eig(A0))/2/pi
A=[zeros(sA,sA) eye(sA,sA);-A0 -A1];
%eig(A(1:sA,1:sA))
[om,z,poles]=damp(A);
nf=om/2/pi;

H=Vscript(H,f*2*pi,0);

%tfplot(f,Vscript(H,f*2*pi,0))
lf=length(f);

Vfit=[Vscript0_1 -Uscript0_1]*AA;
HH=frfgen2(A0,A1,B0,B1,B2,f(lf/2+1:lf));

%tfplot(f(lf/2+1:lf),[HH H(lf/2+1:lf,:)])
%pause
if reviewfit==1
    for iii=1:nmodes
        tfplot(f(lf/2+1:lf),[HH(:,iii) H(lf/2+1:lf,iii)],6)
        legend('Data','Fit')
        figure(gcf)
        clc
        disp(['Reviewing fit ' num2str(iii) ' of ' num2str(nmodes) '. Hit return to continue.'])
        pause
    end
end
%size(H)
u=modefit(f,Horiginal,poles);


%u=modefit(f,Horiginal,poles);

if freqdebug==1, disp('End of Polyref2'),end
return


