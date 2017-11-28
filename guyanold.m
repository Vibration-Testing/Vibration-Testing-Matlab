function [Mr,Kr,T,master,slave]=guyan(m,k,thresh)
%GUYAN Reduce size of second order system of equations by 
% applying Guyan reduction.
% M x'' + K x = 0
% Mr xm'' + Kr xm = 0
% Where x = T*xm, Mr= transpose(T)*M*T, Kr= transpose(T)*K*T
%
%[Mr,Kr,T,master,slave]=GUYAN(M,K,thresh) slaves coordinates i for which 
% thresh>(k(i,i)/m(i,i))/max(k(i,i)/m(i,i)). Master coordinates are returned
% in the vector master, slave coordinate indices are returned in the vector
% slave.
%
%[Mr,Kr,T,master,slave]=GUYAN(M,K,master) slaves coordinates which are not
% listed in the vector master. Master coordinates are returned
% in the vector master, slave coordinate indices are returned in the vector
% slave.
%
%[Kr,T,master,slave]=GUYAN(K,master) performs static condensation, slaving 
% coordinates that are not listed in the vector master. Master coordinates are 
% returned in the vector master, slave coordinate indices are returned in the 
% vector slave.
%
%[Mr,Kr,T,master,slave]=GUYAN(M,K) slaves coordinates i for which 
% 0.1>(k(i,i)/m(i,i))/max(k(i,i)/m(i,i)). Master coordinates are returned
% in the vector master, slave coordinate indices are returned in the vector
% slave.
%
% Reduced coordinate system forces can be obtained by 
% fr=transpose(T)*F
%
% Reduced damping matrix can be obtained using Cr=transpose(T)*C*T.
%
% If mode shapes are obtained for the reduced system, full system mode shapes 
% are phi=T*phi_r
%

%
% Copyright Joseph C. Slater, 6/19/2000
if size(k,1)==1 or size(k,2)==1
	statcond=1;
	thresh=k;
	k=m;
	m=eye(size(k));
end


sprse=1;
if ~issparse(m)
	m=sparse(m);
	k=sparse(k);
	sprse=0;
end
if nargin==2
	thresh=.1;
end

ncoord=length(m);

if length(thresh)==1
	if thresh>=1
		thresh=.1;
		warndlg({'thresh must be less than 1.';'thresh has been set to .1'},'Threshold too hign')
%        disp('move on')
	end
	dm=diag(m);
	dk=diag(k);
	rat=dm./dk;
	%pause
	%plot(sort(1./rat))
	%pause
	mr=rat./max(rat);
	%plot(sort(mr))
	%pause
	mth=min(rat)/max(rat);
	master=(mr)>thresh;
	numofmasters=sum(master);
	[rat,i]=sort(mr);
	slave=i(1:ncoord-numofmasters);
	master=i(ncoord-numofmasters+1:ncoord)';
	lmaster=length(master);
else
	master=thresh;
	i=1:ncoord;
	lmaster=length(master);
	i(master)=zeros(lmaster,1);
	i=sort(i);
	slave=i(lmaster+1:ncoord);
end
if lmaster==ncoord
	figure(gcf)
	plot((1:ncoord)',rat./max(rat),[1 ncoord],[thresh thresh])
	axis([1,ncoord, 0,1])
	legend('Normalized Ratios','Cutoff for Slave Coordinates',0)
	grid on
	xlabel('Coordinate number')
	ylabel('Normalized diagonal ratio.')
	h=warndlg({['thresh must be greater than ' num2str(mth) ' for']...
	;'this problem. Run aborted.'},'Threshold too hign');
	t=(1:.1:70)';
	noise=sin(10*t)+sin(3*t);
	quiet=0*t;
	sound([noise;quiet;noise;quiet;noise])
	slave=0;master=0;Mr=0;Kr=0;T=0;
	break
end
kmm=k(master,master);
ksm=k(slave,master);
kss=k(slave,slave);
kms=ksm';
T=zeros(ncoord,lmaster);
T=sparse(T);
T(master,1:lmaster)=eye(lmaster,lmaster);
T(slave,1:lmaster)=-kss\ksm;
if statcond~=1
	Mr=T'*m*T;
end
Kr=T'*k*T;
if sprse==0
	Mr=full(Mr);
	Kr=full(Kr);
	T=full(T);
end
if statcond==1
	Mr=Kr;
	Kr=T;
	T=master;
	master=slave;
end
