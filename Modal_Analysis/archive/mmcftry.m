% Test case with repeated modes.
clear all
  M=eye(3);
  K=[3 -1 -1;-1 3 -1;-1 -1 3];
  C=.01*K;
  freqs=linspace(0,.5,800)';%+.0001;
  n=length(freqs);
  noise=0.00;%1;
  H=[];
  for el=1:3
	  for em=1:3
		  [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,el,em,freqs);
		  %figure(1)
		  el,em
		  R(:,el,em)=Recep+0*noise*(randn(n,1)+randn(n,1)*i);
		  tfplot(Freq,R(:,el,em));figure(1);pause
		  H=[H Recep];
	  end
  end
  %pause
  [z,nf,a]=mmcf(Freq,R)
damp(state(M,C,K))
