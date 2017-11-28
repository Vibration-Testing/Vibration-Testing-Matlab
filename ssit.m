function [lambda,phi]=ssit(M,K,p,tol,q)
%  SSIT Subspace iteration
%
%  [lambda,phi] = SSIT(M,K,p,tol,q) produces matrices lambda and phi
%  corresponding to the first p eigenvalues and eigenvectors of 
%  the eigenvalue problem M*phi=K*lambda*phi. The variable tol defines 
%  the convergence tolerance such that lambda(i) has converged to tol 
%  places for 1<i<p. tol=6 by default. Note that the tolerance of the 
%  eigenvectors is accurate to tol/2 places. The variable q 
%  defines the number of basis vectors to use and is optional 
%  with the default value of min(2*p,p+8).

% Copyright Joseph C. Slater, May 25, 1998
% Algorithm from Bathe and Wilson, 1976
% Incorporated into Vibration Toolbox 9/23/98

if nargin==3, 
	q=min(2*p,p+8);
	tol=6;
end
if nargin==4
	q=min(2*p,p+8);
end
if p>6
	n=6;
else
	n=1;
end

shift=mean(diag(K))/mean(diag(M))/10^6;
K=K+M*shift;
% Sparse setup
M=sparse(M);
K=sparse(K);
n=length(M);
R=chol(M);%whos
K=R'\K/R;
M=diag(sparse(ones(n,1)));
% Set up initial vectors
MX=zeros(n,q);MX=sparse(MX);
[y,i]=sort(-1./diag(K));
MX(i(1:q-1),1:q-1)=eye(q-1,q-1);
MX(:,q)=diag(M);
X=K\MX;

% We'll do this a lot, so let's get it right now.
kim=K\M;


% We need something to compare to the first time.
% The negative assures convergence won't happen the first time.
%pause
lambdaold=-eye(p,p);

for fa=1:10000
	kp=X'*K*X;%disp('kp')
	mp=X'*M*X;%disp('mp'),
	[phip,lambda]=eig(full(kp),full(mp));%disp('eig')
	[y,i]=sort(diag(lambda));
	phip=sparse(phip(1:q,i));
	lambda=sparse(diag(y));%whos
	disp(['Iteration ' num2str(fa) '. Current 1st natural frequency ' num2str(sqrt(lambda(1,1)-shift)/2/pi) ' Hz.' ])
	X=X*phip;
%   Normalization that is not really needed except for final result
	nx=diag(1./sqrt(sum(X.^2,1)));
	X=X*nx;
	lambda
	% Check to see if first p eigenvalues have converged.
	abs(diag(lambdaold(n:p,n:p)-lambda(n:p,n:p))/diag(lambdaold(n:p,n:p)))
	if (abs(diag(lambdaold(n:p,n:p)-lambda(n:p,n:p))/diag(lambdaold(n:p,n:p))))...
		<10^(-tol)
		break
	end
% 	pause
	X=kim*X;
	lambdaold=lambda;

end
nx=diag(1./sqrt(sum(X.^2,1)));
phi=X*nx;
j;
phi=R\phi;
lambda=lambda-shift*eye(size(lambda));
