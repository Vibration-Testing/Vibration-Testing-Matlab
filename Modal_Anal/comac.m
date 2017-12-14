function comac=comac(u1,u2)
% MAC=MAC(u1,u2)
% Calculate the Coordinate Modal Assurance Criterion for 2 sets of modes
%
% See Harris pp. 21.72

if size(u1)~=size(u2)
	disp('Either the number of modes doesn''t match or the length of the model vectors doesn''t match. Cannot compute.')
	return

else
	mc=mac(u1,u2,1); 
	[y,i]=max(abs(mc)+abs(mc'));
	r=1:size(u1,2);
	comac=zeros(size(u1,1),1);
	%Resort modes to correlate. 
	u2p=u2(:,i);
	mc=mc(:,i);
	u2p=u2p*diag(sign(diag(mc)));
	%[u1 u2p];
	for p=1:size(u1,1)
		comac(p)=sum(u1(p,r).*u2p(p,r))^2/(sum(u1(p,r).*u1(p,r))*sum(u2p(p,r).*u2p(p,r)));
	end
end
