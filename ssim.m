function [xout,xdotout,xddotout]=ssim(M,D,K,f,t,x0,v0)
%SSIM Linear second order system simulation.
% [X,V,A] = SSIM(M,D,K,f,t,x0,v0) simulates the system described by
%      2
%     d x      dx   
%   M --   + D -- + K x = f
%       2      dt
%     dt 
%
% to the force f. The matrix f must have as many columns as
% there are degrees of freedom. Each row "i" of f corresponds to
% the forcing function vector evaluated at t(i). The length of f
% must be the same as the length of t. The vectors x0 and v0
% represent the initial conditions, if they exist.
%
% SSIM(M,D,K,f,t) will plot the response of the system.
%
% Example:
% m=eye(2);                      % mass matrix                                            
% k=[2 -1;-1 1];                 % stiffness matrix                                   
% d=.05*k;                       % damping matrix                                          
% n=4048;dt=.2;                  % 
% t=(0:n-1)'*dt;                 % time vector
% f=[sin(10*t) zeros(size(t))];  %                                
% ssim(m,d,k,f,t)
%
% See also LSIM

% Copyright (c) 1995 by Joseph C. Slater

sm=size(M);

if nargin==5
  x0=zeros(sm(1),1);
  v0=x0;
end

A=[zeros(sm) eye(sm);-M\K -M\D];

B=eye(sm*2);

dt=t(2)-t(1);

Ad=expm(A*dt);
Bd=A\(Ad-eye(size(Ad)))*B;

sf=size(f);
z=ltitr(Ad,Bd,[zeros(size(f)),f],[x0;v0]);

X=z(:,1:sf(2));

V=z(:,sf(2)+1:sf(2)*2);

A=M\(f'-D*V'-K*X');

if nargout==0
  
  plot(t,X)

  title('Displacements versus Time')
  xlabel('Time')
  ylabel('Displacements')
  grid
  zoom on

  pause
  plot(t,V)

  title('Velocities versus Time')
  xlabel('Time')
  ylabel('Velocities')
  grid
  zoom on
  pause
  plot(t,A)

  title('Accelerations versus Time')
  xlabel('Time')
  ylabel('Accelerations')
  grid
  zoom on
  return
end
xout=X;
xdotout=V;
xddotout=A;


