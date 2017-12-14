zeta=.0001;
om=1;
amp=10;
i=sqrt(-1);
w=(0:.01:2)';
f=w/2/pi;
R=amp./(-w.^2+2*i*zeta*w*om+om^2);
[z,nf,a]=sdofcf(f,R)
