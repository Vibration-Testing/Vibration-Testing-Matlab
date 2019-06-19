% Combine multiple signals into just one. Works great.

  M=eye(2);
  K=[2 -1;-1 2];
  C=.01*K;
  [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,2,linspace(0,.5,1024));
  figure(1)
  n=1*250;
  m=.5;
  f=Freq((1:n)+450*m);
  R2=Recep((1:n)+450*m);
  R2=R2+.1*randn(n,1)+.1*randn(n,1)*i;% Poorly Simulated Noise
  [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,1,linspace(0,.5,1024));
  R1=Recep((1:n)+450*m);
  R1=R1+.1*randn(n,1)+.1*randn(n,1)*i;% Poorly Simulated Noise
  [z,nf,a]=mdofcf(f,[R1 R2])

