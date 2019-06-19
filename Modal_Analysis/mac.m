function ModeAC=mac(u1,u2,opt)
% MAC=MAC(Phi1,Phi2)
% Calculate the Modal Assurance Criterion for 2 sets of modes
%
% Output is a matrix such that the n, m element correponds to the MAC value
% between the nth mode of Phi1 and the mth mode of Phi2. 

ModeAC=[];

if size(u1)~=size(u2)
	
	disp('Either the number of modes doesn''t match or the length of the model vectors doesn''t match. Cannot compute.')
ModeAC=0;
return
else
	if nargin==2
	for i=1:size(u1,2)
		for j=1:size(u1,2)
            %u1(:,i)'*u2(:,j);
			ModeAC(i,j)=(u1(:,i)'*u2(:,j))^2/((u1(:,i)'*u1(:,i))*(u2(:,j)'*u2(:,j)));%
		end
	end
	else
		for i=1:size(u1,2)
		for j=1:size(u1,2)
			ModeAC(i,j)=(u1(:,i)'*u2(:,j))/sqrt(((u1(:,i)'*u1(:,i))*(u2(:,j)'*u2(:,j))));%
		end
		end
	end
end
ModeAC=ModeAC';