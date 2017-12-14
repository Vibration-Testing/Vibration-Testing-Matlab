%complex mode indicator test script
clear 

k=[2 -1;-1 2];
m=eye(size(k));
c=.01*k;

for i=1:2
	for j=1:2
	[w,R,M,I]=vtb7_5(m,c,k,i,j,linspace(0,.5,1024));
	H(i,j,:)=R;
	end
end

%cmif(w,H)
clear 
%pause
k=[3 -1 -1;-1 3 -1;-1 -1 3];
m=eye(size(k));
sqrt(eig(m\k))
c=.01*k;
hold off
for i=1:3
	for j=1:3
	[f,R,M,I]=vtb7_5(m,c,k,i,j,linspace(0,2/pi,1024));
	H(i,j,:)=R+.0*randn(size(R));
	%tfplot(f,R)
	%pause
	end
end
size(H);
cmif(f*2*pi,H)

