
%VTDEM Demo for Frequency Domain Section of the Vibration Testing Toolbox
% Run by typing sshow('vtdem')
%

% Copyright Joseph C. Slater, 1995

%
numdat=5; %number of runs for averaging
clear y u
%if ~exist('SlideShowGUIFlag'),sshow('vtdem'),return;end
if ~exist('SlideShowGUIFlag'), figNumber=0;end

load vtdem1
subplot(1,1,1)
imhan=image(x);
colormap(map)
axis('image')
axis('off')
cfun=get(gcf,'units');
set(gcf,'units','pixels')
cfpos=get(gcf,'position');
ScreenUnits=get(0,'units');
set(0,'units','pixels');
Screensize=get(0,'screensize');
height=344;width=428;
%cfpos(1:2)=[cfpos(1) cfpos(2)+cfpos(4)-height];
%cfpos(3:4)=[width height];
cfpos(1)=(Screensize(3)-width)/2;
cfpos(2)=(Screensize(4)*5/6-height);
cfpos(3)=width;
cfpos(4)=height;
set(gcf,'position',cfpos)
set(gcf,'units',cfun,'resize','off')

if ssinit(figNumber)
  str=...
   ['Let''s consider taking a simple two degree of freedom     '
    'spring-mass-damper system connected to ground at the     '
    'left end. Each mass is 1 kg, the springs have a stiffness'
    'of 1 N/m, and dashpots run parallel to the springs       '
    'with damping of .05 N-s/m.                               '
    '>>m=eye(2);                                              '
    '>>k=[2 -1;-1 1];                                         '
    '>>d=.05*k;                                               '];
  ssdisp(figNumber,str);
  if figNumber,return;end
end
m=eye(2);
k=[2 -1;-1 1];
d=.05*k;

str=...
 ['The transfer function between a force at degree    '
  'of freedom 1 and the motion at degree of freedom 1 '
  'can be observed using TF. We''ll plot only the      '
  'compliance transfer function here, but the mobility'
  'and inertance transfer functions can also be easily'
  'obtained using TF.                                 '
  '>>sostf(M,D,K,1,1,[0 1])                           '];
  
ssdisp(figNumber,str);

[Freq,Recep,Mobil,Inert]=sostf(m,d,k,1,1,[0 1]);
tfplot(Freq,Recep)
%subplot(2,1,1)
title('Receptance')
%delete(imhan);
set(gcf,'resize','on')
chil=get(gcf,'children');

if sspause(figNumber),return;end;

dt=.2;
n=4048;
t=(0:n-1)'*dt;
u=randn(n,1);

for i=1:length(chil)
  if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
    delete(chil(i))
	i=length(chil);
	drawnow
  end
end

str=...
 ['To simulate the system under random excitation, we first  '
  'generate a time vector, then use the function RANDN to    '
  'generate a random signal.                                 '
  '>>dt=.2;n=4048;                                           '
  '>>t=(0:n-1)''*dt;                                          '
  '>>u=randn(n,1);                                           '];

ssdisp(figNumber,str);
subplot(1,1,1)
curaxe=gca;
plot(t,u),
grid on,
xlabel('Time'),
ylabel('Forcing Function'),
title('Random Force')
%delete(curaxe)

if sspause(figNumber),return;end;

str=...
 ['It can be seen that this forcing function has a good   '
  'distribution as well (zero mean, nearly Gausian).      '
  '>>hist(u);                                             '
  '>>grid                                                 '];

ssdisp(figNumber,str);
curaxe=gca;
delete(curaxe)
hist(u),grid
axis([-4 4 0 1500])
if sspause(figNumber),return;end;

str=...
 ['This forcing function has a near constant auto      '
  '(power) spectral density too...                     '
  '>>asd(u,t)                                          '];
  
ssdisp(figNumber,str);
asd(u,t)
chil=get(gcf,'children');
for i=1:length(chil)
  if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
    delete(chil(i))
	i=length(chil);
	drawnow
	break
  end
end

if sspause(figNumber),return;end;

str=...
 ['For a random signal, the auto correlation should be '
  'zero everywhere except tau=0.                       '
  '>>crcor(u,u,dt)                                     '];
  
ssdisp(figNumber,str);
crcor(u,u,dt) 
chil=get(gcf,'children');
for i=1:length(chil)
  if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
    delete(chil(i))
	i=length(chil);
	drawnow
	break
  end
end

for i=2:numdat
  u1=randn(length(t),1);
  u=[u u1];
end

if sspause(figNumber),return;end;

str=...
 ['Running this process repeatedly, we can obtain a set of'
  '5  forcing functions to excite the system with.        '
  '>>for i=2:5                                            '
  '>>  u1=randn(length(t),1);                             '
  '>>  u=[u u1];                                          '
  '>>end                                                  '
  '>>plot(t,u)                                            '];

ssdisp(figNumber,str);
plot(t,u)
grid on,
xlabel('Time'),
ylabel('Forcing Function'),
title('Random Forces')
pause(.0001)


if sspause(figNumber),return;end;

str=...
 ['The cross correlation between each of these signals    '
  'should be zero everywhere if they are truly random.    '
  'This can be observed using the function CRCOR on any   '
  'pair of signals.                                       '
  '>>crcor(u(:,1),u(:,2),dt)                              '
  '>>axis([-1000 1000 -1000 5000]                         '];

ssdisp(figNumber,str);

[Tau,Pxy]=crcor(u(:,1),u(:,2),dt);
plot(Tau,Pxy)
title('Linear Cross Correlation')
xlabel('Tau')
ylabel('Linear Cross Correlation')
grid
axis([-1000 1000 -1000 5000])
drawnow

zoom on
xtot=[];
for i=1:numdat
  [x,v,a]=ssim(m,d,k,[u(:,i) zeros(length(u(:,1)),1)],t);
  xtot=[xtot x(:,1)];
end


if sspause(figNumber),return;end;

str=...
 ['The responses of the first mass to these excitations        '
  'can be obtained using the function SSIM.                    '
  '>>for i=1:5                                                 '
  '>>  [x,v,a]=ssim(m,d,k,[u(:,i) zeros(length(u(:,1)),1)],t); '
  '>>  y=[y x(:,1)];                                           '
  '>>end                                                       '];

ssdisp(figNumber,str);
y=xtot;
plot(t,y)
grid on
xlabel('Time')
ylabel('Response')
title('Response of a 2DOF system to a random excitation')
[Freq,Txf]=tfest(y,u,t);

if sspause(figNumber),return;end;


str=...
 ['Unlike the forcing functions, the auto correlation of'
  'the response is not zero everywhere other than tau=0.'
  'The average auto correlation shows that the signals  '
  'are not completely random.                           '
  '>>crcor(y,y,t)                                       '];

crcor(y,y,t)
ssdisp(figNumber,str);
chil=get(gcf,'children');
for i=1:length(chil)
  if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
    delete(chil(i))
	i=length(chil);
	drawnow
	break
  end
end

if sspause(figNumber),return;end;

str=...
 ['Looking at the average power spectral density using CRSD      '
  'demonstrates that the response does indeed contain significant'
  'energy at two distinct frequencies.                           '
  '>>crsd(y,y,t)                                                 '];

ssdisp(figNumber,str);
crsd(y,y,t)
title('Auto Spectrum Density')
chil=get(gcf,'children');
for i=1:length(chil)
  if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
    delete(chil(i))
	i=length(chil);
	drawnow
	break
  end
end
if sspause(figNumber),return;end;

str=...
 ['An estimate of the transfer function can be obtained using the'
  'function TFEST. You can zoom in on the plot by clicking on a  '
  'point in the plot and zoom out by using the right mouse button'
  '(shift click on the Mac). You can also drag a box around an   '
  'area you''d like to zoom into. Double clicking brings you back '
  'to the full view. TFEST will automatically average the data   '
  'transfer for you by default. You can specify no averaging and '
  'modify the default zero padding as well.                      '
  '>>tfest(y,u,t);                                               '];

ssdisp(figNumber,str);
tfplot(Freq,Txf)
axis([0 2.5 -240 240])
%pause
% chil=get(gcf,'children');
% for i=1:length(chil)
%   if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
%     delete(chil(i))
% 	i=length(chil);
% 	drawnow
% 	break
%   end
% end
if sspause(figNumber),return;end;

str=...
 ['The actual estimate can be obtained in vector form by'
  'using left hand side arguments. The results may then '
  'be plotted using TFPLOT. TFPLOT gives you much more  '
  'control over the way the transfer function estimate  '
  'is displayed. For instance, let''s look at the        '
  'frequency range near the peaks...                    '
  '>>[Freq,Txf]=tfest(y,u,t);                           '
  '>>tfplot(Freq,Txf,0,.4)                              '];

ssdisp(figNumber,str);
tfplot(Freq,Txf,0,.6)
% chil=get(gcf,'children');
% for i=1:length(chil)
%   if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
%     delete(chil(i))
% 	i=length(chil);
% 	drawnow
% 	break
%   end
% end
if sspause(figNumber),return;end;

str=...
 ['or perhaps the Nyquist plot...                      '
  '>>tfplot(Freq,Txf,0,.6,5)                           '];

ssdisp(figNumber,str);
tfplot(Freq,Txf,0,.6,5)
% chil=get(gcf,'children');
% for i=1:length(chil)
%   if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
%     delete(chil(i))
% 	i=length(chil);
% 	drawnow
% 	break
%   end
% end
if sspause(figNumber),return;end;
tic
tfplot(Freq,Txf,0,.6)
% chil=get(gcf,'children');
% for i=1:length(chil)
%   if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
%     delete(chil(i))
% 	i=length(chil);
% 	drawnow
% 	break
%   end
% end

str=...
 ['These results do not look very good, though.  This is'
  'likely the results of leakage (using incomplete data '
  'cycles). We can improve the results by applying a    '
  'window of some form, for example a Hanning (von      '
  'Hanning) window.                                     '
  '>>plot(t,hannwin(u))                                 '
  '>>plot(t,hannwin(y))                                 '];

ssdisp(figNumber,str);
drawnow
uh=hannwin(u);
yh=hannwin(y);
%pause(5-toc)
subplot(2,1,1)
plot(t,uh)
grid on
xlabel('Time')
%ylabel('Windowed Forcing Function')
title('Windowed Forcing Function')
subplot(2,1,2)
plot(t,yh)
grid on
xlabel('Time')
%ylabel('Windowed Response')
title('Windowed Response')
drawnow
[Freq,Txf]=tfest(hannwin(y),hannwin(u),t);
if sspause(figNumber),return;end;

str=...
 ['After applying a von Hanning window to the forcing '
  'function and to the response, the transfer function'
  'estimate looks much better.                        '
  '>>[Freq,Txf]=tfest(hannwin(y),hannwin(u),t);       '
  '>>tfplot(Freq,Txf,0,.6)                            '];

ssdisp(figNumber,str);

tfplot(Freq,Txf,0,.6)
% chil=get(gcf,'children');
% for i=1:length(chil)
%   if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
%     delete(chil(i))
% 	i=length(chil);
% 	drawnow
% 	break
%   end
% end

if sspause(figNumber),return;end;

str=...
 ['The quality of the transfer function can be    '
  'verified by plotting the coherance function    '
  'using COH.                                     '
  '>>coh(y,u,t)                                   '];

ssdisp(figNumber,str);

subplot(1,1,1)
coh(hannwin(y),hannwin(u),t) 
axis([0 .6 0 1])

% chil=get(gcf,'children');
% for i=1:length(chil)
%   if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
%     delete(chil(i))
% 	i=length(chil);
% 	drawnow
% 	break
%   end
% end

if sspause(figNumber),return;end;

str=...
 ['These results can be compared to the "true" transfer          '
  'obtained earlier by using TFPLOT to plot them over each other.'
  '>>tfplot(Freq,[Txf, Recep],0,.5)                              '];

ssdisp(figNumber,str);

[Freq,Recep,Mobil,Inert]=sostf(m,d,k,1,1,Freq);

tfplot(Freq,[Txf, Recep],0,.5)
% chil=get(gcf,'children');
% for i=1:length(chil)
%   if strcmp(get(chil(i),'type'),'axes') & chil(i)~=gca
%     delete(chil(i))
% 	i=length(chil);
% 	drawnow
% 	break
%   end
% end

