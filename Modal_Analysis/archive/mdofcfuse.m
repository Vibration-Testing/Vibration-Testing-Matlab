  clf
  M=eye(2)*2;
  K=[2 -1;-1 2];
  C=.01*K;
  [Freq,Recep1,Mobil1,Inert1]=vtb7_5(M,C,K,1,1,linspace(0,.5,1024)); %Frequency response function between node 1 and itself
  [Freq,Recep2,Mobil2,Inert2]=vtb7_5(M,C,K,1,2,linspace(0,.5,1024)); %Frequency response function between node 1 and node 2
  figure(1)
  plot(Freq,20*log10(abs(Recep1)))
  Recep=[Recep1 Recep2];
  Recep=Recep+.1*randn(size(Recep))+i*.1*randn(size(Recep));
  [z,nf,u]=mdofcf(Freq,Recep,.1,.12);
  z(1,1)=z;
  lambda(1,1)=(nf*2*pi)^2;
  
  S(:,1)=real(u);%/abs(sqrt(u(1)));
  [z,nf,u]=mdofcf(Freq,Recep,.16,.25);
  z(2,2)=z;
  lambda(2,2)=(nf*2*pi)^2;
  S(:,2)=real(u);%/abs(sqrt(u(1)));
  dampingratios=diag(z)
  naturalfrequencies=sqrt(diag(lambda))/2/pi
  M=M
  Midentified=S'\eye(2)/S%Make Mass matrix
  K=K
  Kidentified=S'\lambda/S%Make Stiffness Matrix
  C=C
  Cidentified=S'\(2*z*sqrt(lambda))/S%Make damping matrix
  
