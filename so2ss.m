function [a,b,c,d]=so2ss(M, D, K, Bt, Cd, Cv, Ca)
% SO2SS [a,b,c,d]=SO2SS(M, C, K, Bt, Cd, Cv, Ca) returns the state system
% [xdot ]   [   0        I  ][ x  ]    [   0   ]
% [     ] = [               ][    ] +  [       ] [u]
% [xddot]   [-M^-1*K -M^-1*C][xdot]    [M^-1*Bs]
%
% y = C z + D u
% 
% where z = [x' xdot']'
%
% given the second order form equations
% Mxddot+Cxdot+x=Bs u
% y = Cd x + Cv xdot + Ca xddot
  
% Copyright Joseph C. Slater, 2003
% Updates Dec 11, 2017 to include Cd, Cv, and Ca.

l=size(M,1);

if nargin==4
    Cd = eye(l)
    Ca = zeros(l)
    Cv = Ca
if nargin<4
  Bt=eye(l)  
if nargin<3
  K=D;D=0*K;
end

a=[zeros(l) eye(l);-M\K -M\D];
b=[zeros(l,size(Bt,2)); M\Bt];
c=[Cd-Ca*M\K, Cv-Ca*M\D]
d=Ca*M\Bt;

