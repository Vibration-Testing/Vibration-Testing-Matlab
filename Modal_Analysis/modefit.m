function [u]=modefit(f,H,poles)
% f is in fricking Hz!
%disp('Extracting modes')
u=[];
if size(f,1)==1
    f=f';
end
size(H);
if size(poles,2)==1
    poles=poles.';
end

p=size(H,1)/length(f);
if length(f)~=size(H,1)
    if p~=floor(p)
        disp(['You cannot have ' num2str(size(H,2)/length(f)) ' outputs'])
        return
    end    
%    H=Vscript(H,f*2*pi,0);
%    size(H);
end
ln=length(f)/2;
%f=[-f;f];
%H=[conj(H);H];


modalh=zeros(size(H,1),length(poles)+2);

poles;

[P,OM]=meshgrid(poles,2*pi*f*1i);
%abs(poles)/2/pi
%poles
%size(H)
modalh(1:p:size(H,1),:)=[-1./(2*pi*f).^2 ones(length(f),1) 1./(OM-P)];
%disp('size modalh')
%size(modalh)
%tfplot(f,[modalh conj(modalh)])
for ii=2:p
    modalh(ii:p:size(H,1),:)=modalh(1:p:size(H,1),:);
end
modalh;

BB=modalh\H;
lp=length(poles);
%B=transpose(BB(3:2:size(BB,1),:))
%diag(1./(poles(1:2:lp)-conj(poles(1:2:lp))))
u=transpose(BB(3:2:size(BB,1),:))*diag(1./(poles(1:2:lp)-conj(poles(1:2:lp))));

size(f);
ln;
ff=f((ln+1):2*ln,:);
HH=H(ln+1:2*ln,:);
tfplot(ff,[modalh(ln+1:ln*2,:)*BB HH])
tfplot(ff,(H))
%plot(f(ln+1:2*ln,:),20*log10(abs(H(ln+1:2*ln,:))))
%plot(f(ln+1:2*ln,:),f(ln+1:2*ln,:))
max(f(ln+1:2*ln,:));
min(f(ln+1:2*ln,:));
