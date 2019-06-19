M=eye(2);K=[2 -1;-1 2];C=.01*K;
%[Freq,Recep,Mobil,Inert]=sostf(M,C,K,1,1,linspace(0,.5,1024));
[Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,1,linspace(0,.5,1024));

figure(1)
%tfplot(Freq,Recep)
n=250;
f2=Freq((1:n)+450);
R2=Recep((1:n)+450);
R2=R2+.1*randn(n,1)+.1*randn(n,1)*i;
[z,w,a,com]=vtb7_4(f2,R2)
