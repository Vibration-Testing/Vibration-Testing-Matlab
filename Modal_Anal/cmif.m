function CMIF(w,H)
% CMIF(w,H) calculates and plots the complex mode indicator function
%  For the FRFs, H, rows and columns reference degrees of freedom, and
%  the third dimension is wrt omega.

% Copyright Joseph C. Slater, 2001
% Version 1.0

s=size(H);
if min(s(1:2))==1
	warndlg('Without MIMO data, this function means little','Huh?')
end
l=s(3);
w=w(:);
if s(2)>s(1)
	a=H(:,:,1)*H(:,:,1)';
	[v(:,:,1),lam(:,:,1)]=eig(real(a));
	for i=1:l
	  H(:,:,i);
	  a=H(:,:,i)*H(:,:,i)';
	  a=v(:,:,i-1)\a*v(:,:,i-1)
	  pause
	  [v(:,:,i),lam(:,:,i)]=eig(real(a));
	  v(:,:,i);
  end
else
	a=H(:,:,1)'*H(:,:,1);
	[v(:,:,1),lam(:,:,1)]=eig(real(a));
	for i=2:l
	  a=real(H(:,:,i)'*H(:,:,i));
	%a=v(:,:,i-1)\real(a)*v(:,:,i-1);
	%pause
	  [v(:,:,i),lam(:,:,i)]=eig(a);
	  v(:,:,i)\real(a)*v(:,:,i);
	  v(:,:,i);
  end
end
size(lam);
s;
for i=1:min(s(1:2))
	size(lam(i,i,:));
	d=lam(i,i,:);
	size(d);
	d=d(:);
	size(d);
%	disp('size d');
	y(:,i)=d;
end
y;
size(y)
%y(20,:)
y=-(sort(-y'))';
%semilogy(w,y(:,1),'>',w,y(:,2),'+',w,y(:,3),'o')
%axes

%set(0,...%'DefaultAxesColorOrder',[0 0 1],...
%             'DefaultAxesLineStyleOrder','o|x|+|*|S|D|v|^|<|>|p|h')
%set(0,'DefaultAxesLineStyleOrder','o|x|+|*|S|D|v|^|<|>|p|h')
%'yo|mx|c+|r*|gS|bD|wv|k^|k<|k>|kp|kh')
ls=['bo';'gx';'r+';'c*';'mS';'kD';'bv';'g^';'r<';'c>';'mp';'kh'];
semilogy(w,y(:,1),ls(1,:))
hold on
for i=1:size(y,2)
semilogy(w,y(:,i),ls(i,:))
end
hold off
grid on
zoom on
