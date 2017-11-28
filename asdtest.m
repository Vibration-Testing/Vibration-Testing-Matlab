t=0:1/1024:1-1/1024;
a=sin(t*2*pi*20);
t2=0:1/256:1-1/256;
a2=sin(t2*2*pi*20);
[f,p]=asd(a',t);
[f2,p2]=asd(a2',t2);
semilogy(f,p,f2,p2)
