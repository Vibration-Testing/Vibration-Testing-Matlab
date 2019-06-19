global w2 R2 XoF
M=eye(2);K=[2 -1;-1 2];C=.01*K;
[Freq,Recep,Mobil,Inert]=sostf(M,C,K,1,1,linspace(0,.5,1024));
figure(1)
tfplot(Freq,Recep)
n=350*2
f2=Freq(1:n);
R2=Recep(1:n);
R2=R2+.1*randn(n,1)+.1*randn(n,1)*i;
figure(2)
tfplot(f2,R2)
w2=2*pi*f2;
format long

x0=[.5;.005;1;-1;.01;1.6;0;0;0];x=x0;
%x0=[.5;.005;1;0;0;0];x=x0;
%plot(f2,20*log10(abs(x(1)./(-w2.^2+2*x(2)*w2*i*x(3)+x(3)^2)+x(4)+x(5)*i*w2-x(6)*w2.^2)),'g',f2,20*log10(abs(R2)))
x0

clear i
for j=1:100

x=fmins('cferror',x0)

figure(3)
plot(f2,20*log10(abs(XoF)),'g',f2,20*log10(abs(R2)))
zoom on
                   % x(1)./(-w2.^2+2*x(2)*w2*i*x(3)+x(3)^2)+x(4)+x(5)*i*w2-x(6)*w2.^2;
grid on
%hold on
%figure(3)
%plot(f2,20*log10(abs(R2)))
%hold off
x0=x;
pause
end
