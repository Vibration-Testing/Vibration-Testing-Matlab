function femgui(mode)
% femgui  [node,ncon,zero,force]=femgui
% Makes input file for VTB8_2
% I suggest you print the file vtb8read.

%	Joseph C. Slater, 6-17-90
%	Copyright (c) 1990-94 by Joseph C. Slater

%       21 Dec 00 -- Major revision of GUI.
%       11 Nov 94 -- fixed graphics to work with Matlab 4.0

% What we need to do is instead of having matrices of locations...
% is have matrices of handles for objects whose properties are available
% in the userdata. This way they can each be selected and edited one at a time. 
% Only when writing the data file or reading it to we actually transfer between the two. 
% We could actually set a flag that tells the callback for a node whther to 
% choose that node for deleting, editing, connecting with an element, adding bc, or other...
% In the file menu, we can have an open for edit, or a regular open command.
% The first would simply cause a text edit of the datafile. If no data existed, one would have
% to be written first. Otherwise, one could be selected. However, this option would only be available
% after selection of the new or open commands.
% Callback functions of all object must be turned on and off depending on what objects are currently 
% being edited.


if nargin==0
	clear
	global node
	h=findobj('tag','femgui');
	if isempty(h)
		fn=figure('Name', 'Engineering Vibration Toolbox FEM', ...
	    'NumberTitle', 'off', ...
	    'MenuBar', 'none', ...
	    'Units', 'normalized', ...
	    'PaperOrientation', 'Landscape', ...
	    'PaperUnits', 'Normalized', ...
		'tag','femgui',...
	    'ShareColors', 'No');
	drawnow
	
	file_opt = uimenu(fn,'Label', 'File ');
    label = str2mat(...%label,...
    'New',...
    'Open',...
    'Save',...
	'Analyze',...
	'Quit');
	file_new_m = uimenu(file_opt, ...
    'Label', deblank(label(1,:)), ...
    'Accelerator','n',...
    'CallBack', 'femgui(''new'');');
	file_open_m = uimenu(file_opt, ...
    'Label', deblank(label(2,:)), ...
    'Accelerator','o',...
    'CallBack', 'femgui(''open'');');
	file_save_m = uimenu(file_opt, ...
    'Label', deblank(label(3,:)), ...
    'Accelerator','s',...
	'enable','off',...
    'CallBack', 'femgui(''save'');');
	file_analyze_m = uimenu(file_opt, ...
    'Label', deblank(label(4,:)), ...
    'Accelerator','a',...
	'enable','off',...
    'CallBack', 'femgui(''analyze'');');
	file_quit_m = uimenu(file_opt, ...
    'Label', deblank(label(5,:)), ...
    'Accelerator','w',...
    'CallBack', 'femgui(''quit'');');
	
	file_opt = uimenu(fn,'Label', 'Edit');
    label = str2mat(...%label,...
    'Nodes',...
    'Elements',...
	'Point Masses/Inertias',...
	'Boundary Conditions');
	file_open_m = uimenu(file_opt, ...
    'Label', deblank(label(1,:)), ...
	'enable','off',...
    'CallBack', ['global node;femgui(''nodes'');']);
	file_open_m = uimenu(file_opt, ...
    'Label', deblank(label(2,:)), ...
	'enable','off',...
    'CallBack', 'femgui(''elements'');');
	%get(file_open_m)
	file_open_m = uimenu(file_opt, ...
    'Label', deblank(label(3,:)), ...
	'enable','off',...
    'CallBack', 'femgui(''masses'');');
	file_open_m = uimenu(file_opt, ...
    'Label', deblank(label(4,:)), ...
	'enable','off',...
    'CallBack', 'femgui(''bcs'');');
	%openvar
	%disp('lkhlhl')
	mode='a';	
	else 
		figure(h)
		fn=h;
	end

else
	userdata=get(gcf,'Userdata');
	if ~isempty(userdata)
		global node
		node=userdata{1}
		ncon=userdata{2};
		conm=userdata{3};
		zero=userdata{4};
		force=userdata{5};
	end
	disp('Command')
	if strcmp(mode,'nodes')
		gcbo
		disp('nodes')
		size(node)
		openvar('node')
		disp('get rid of openvar and replace with edit')
		pause(.01),drawnow
		h=questdlg('accept changes?','yes','now')
		node(1:10,1:2)
	elseif strcmp(mode,'quit')
		disp('quit')
		delete(gcf)
	elseif strcmp(mode,'nodes')
		disp('Edit Nodes')
		for i=1:1000
		  loc=input('Enter x and y location of node (ie. [x y]), 0 to end, or 1 to edit nodes. ');
		  if loc==0 & length(loc)==1 ,break,end
		  node(i,:)=loc;
		  length=max([max(node(:,1))-min(node(:,1)) max(node(:,2))-min(node(:,2))]);
		  d=node;
		  xlo=min([node(:,1);d(:,1)]);
		  xho=max([node(:,1);d(:,1)]);
		  ylo=min([node(:,2);d(:,2)]);
		  yho=max([node(:,2);d(:,2)]);
		  xsp=xho-xlo;
		  ysp=yho-ylo;
		  if .618*xsp > ysp
		      xh=xho+.1*xsp;
		      xl=xlo-.1*xsp;
			  yh=(yho+ylo)/2+(xh-xl)*.618/2;
			  yl=(yho+ylo)/2-(xh-xl)*.618/2;
			else
		      yh=yho+.1*ysp;
		      yl=ylo-.1*ysp;
		      xh=(xho+xlo)/2+.5*(yh-yl)/.618;
		      xl=(xho+xlo)/2-.5*(yh-yl)/.618;
		  end
		%%%
		  dx=.01*(xh-xl);
		  dy=.035*(yh-yl);
		  plot(node(:,1),node(:,2),'*b')
		  j=1:i;
		  if i < 10
		      istr=[num2str(i) '  '];
			elseif i < 100
			  istr=[num2str(i) ' '];
			else
			  istr=[num2str(i)];
		  end	
		  nnum=[nnum;istr];
		  text(node(j,1)'+dx,node(j,2)'+dy,nnum)
		 % axis([xl xh yl yh])
		  axis('equal'),grid on
		  al=axis;alx=al(2)-al(1);aly=al(4)-al(3);
		  axis([al(1)-alx/10 al(2)+alx/10 al(3)-aly/10 al(4)+aly/10])
		  clc
	  end
  elseif strcmp(mode,'ed')

  elseif strcmp(mode,'el')
	hold on
	clc
	%disp('Do you have a pointing device such as a mouse or trackball? (y/n)')
	%disp('(Arrow keys may be sufficient)')
	%point=input(' ','s');
	point='y';
	%connecting nodes section
	clc
	home
	  disp(' Pick nodes to connect with elements.')
	if point=='y'
	  disp(' (Use pointing device or arrow keys and return)')
	 else
	  disp(' (Enter node numbers one at a time)')
	end
	  pause(1)
	answer2='n';
	for i=1:1000
	%  clc
	%  home
	if i~=1
	  disp(' Pick nodes to connect with elements.')
	end
	  if point=='y'
	    [x1 y1]=ginput(1);
	    dis=(node(:,1)-x1).^2+(node(:,2)-y1).^2;
	    [dsq,nodenum1]=min(dis);
	   else
	    clc,home
	    nodenum1=input(' Enter node number 1: ');
	  end
	  plot(node(nodenum1,1),node(nodenum1,2),'*w')
	  plot(node(nodenum1,1),node(nodenum1,2),'or')
	  if point=='y'
	    [x2 y2]=ginput(1);
	    dis=(node(:,1)-x2).^2+(node(:,2)-y2).^2;
	    [dsq,nodenum2]=min(dis);
	   else
	    clc,home
	    nodenum2=input(' Enter node number 2: ');
	  end
	  plot(node(nodenum2,1),node(nodenum2,2),'*w')
	  plot(node(nodenum2,1),node(nodenum2,2),'or')
	  pause(.04)
	  plot([node(nodenum1,1) node(nodenum2,1)],[node(nodenum1,2) node(nodenum2,2)],'-b')
	  plot(node(nodenum1,1),node(nodenum1,2),'ow')
	  plot(node(nodenum1,1),node(nodenum1,2),'*b')
	  plot(node(nodenum2,1),node(nodenum2,2),'ow')
	  plot(node(nodenum2,1),node(nodenum2,2),'*b')
	  clc
	  home
	
	  if i>1
	  answer2=input(' Same properties as previous element? (y/n) ','s');
	  end
	  if answer2=='n'
	    clc
	    E=input(' Enter the modulus of elasticity of the member. ');
	    G=input(' Enter the shear modulus of the member (zero for EB beam). ');
	    I=input(' Enter the moment of area of the member. ');
	    A=input(' Enter the cross sectional area of the member. ');
	    Rho=input(' Enter the density per unit length of the member. ');
	  end
	  ncon(i,:)=[nodenum1 nodenum2 E A I G Rho];
	  answer=input(' Enter another element? (y/n) ','s');
	  if answer~='y',break,end
  end
elseif strcmp(mode,'am')
		% adding concentrated masses and inertias
	conm=[];
	for i=1:1000
	clc
	home
	  if i==1
	    answer=input(' Add concentrated masses and rotational inertias? (y/n) ','s');
	  end
	  if i>1
	    answer=input(' Add more concentrated masses and rotational inertias? (y/n) ','s');
	  end
	  if answer~='y',break,end
	  disp(' ')
	  disp(' Pick node to add mass/rotational inertia to.')
	  pause(.5)
	  if point=='y'
	    [x1 y1]=ginput(1);
	    dis=(node(:,1)-x1).^2+(node(:,2)-y1).^2;
	    [dsq,nodenum]=min(dis);
	   else
	    nodenum=input(' Enter node number: ');
	  end
	  plot(node(nodenum,1),node(nodenum,2),'*w')
	  plot(node(nodenum,1),node(nodenum,2),'xr')
	  answer=input(' Add mass or rotational inertia?(m,i,n(one)) ','s');
	  massval=input(' Enter magnitude of mass/inertia. ');
	  conm(i,:)=[0 0 0];
	  conm(i,1)=nodenum;
	
	  if answer=='m'
	    conm(i,2)=massval;
	  end
	  if answer=='i'
	    conm(i,3)=massval;
	  end
	
	  plot(node(nodenum,1),node(nodenum,2),'xi')
	  plot(node(nodenum,1),node(nodenum,2),'*b')
  end
elseif strcmp(mode,'abc')
	
	
	% zeroing of displacements
	zero=[];ij=0;
	for i=1:1000
	clc
	home
	  if i==1
	    answer=input(' Add boundary conditions? (y/n) ','s');
	  end
	  if i>1
	    answer=input(' Zero another displacement? (y/n) ','s');
	  end
	  if answer~='y',break,end
	  disp(' ')
	  disp(' Pick node to zero')
	  pause(1.)
	  if point=='y'
	    [x1 y1]=ginput(1);
	    dis=(node(:,1)-x1).^2+(node(:,2)-y1).^2;
	    [dsq,nodenum]=min(dis);
	   else
	    nodenum=input(' Enter node number: ');
	  end
	  plot(node(nodenum,1),node(nodenum,2),'*w')
	  plot(node(nodenum,1),node(nodenum,2),'xr')
	  disp(' ')
	  disp(' ')
	  disp(' Zero which displacement(s)?(x,y,t(heta),n(one))')
	  disp(' (ie xt for x and theta)');
	  answern=input(' ','s');
	  sanswern=size(answern);
	for ii=1:sanswern(2)
	    ij=ij+1;
	    answer=answern(ii);
	    if answer=='x'
	       plot(node(nodenum,1),node(nodenum,2),'xi')
	       plot(node(nodenum,1),node(nodenum,2),'*b')
	       zero(ij,:)=[nodenum 1];
	    end
	    if answer=='y'
	       plot(node(nodenum,1),node(nodenum,2),'xi')
	       plot(node(nodenum,1),node(nodenum,2),'*b')
	       zero(ij,:)=[nodenum 2];
	    end
	    if answer=='t'
	       plot(node(nodenum,1),node(nodenum,2),'xi')
	       plot(node(nodenum,1),node(nodenum,2),'*b')
	       zero(ij,:)=[nodenum 3];
	    end
	    if answer=='n'
	       plot(node(nodenum,1),node(nodenum,2),'xi')
	       plot(node(nodenum,1),node(nodenum,2),'*b')
	    end
	  end
	end
	elseif strcmp(mode,'af')
	force=[];
	% adding loads
	for i=1:1000
	clc
	home
	  answer=input(' Add loads? (y/n) ','s');
	  if answer~='y',break,end
	  disp(' ')
	  disp(' Pick node to load')
	  pause(1.5)
	  if point=='y'
	    [x1 y1]=ginput(1);
	    dis=(node(:,1)-x1).^2+(node(:,2)-y1).^2;
	    [dsq,nodenum]=min(dis);
	   else
	    nodenum=input(' Enter node number: ');
	  end
	  plot(node(nodenum,1),node(nodenum,2),'*w')
	  plot(node(nodenum,1),node(nodenum,2),'xr')
	  answer=input('Load which displacement?(x,y,t(heta),n(one)) ','s');
	
	  if answer=='x'
	    loadval=input(' Enter magnitude of load. ');
	    force(i,:)=[nodenum 1 loadval];
	  end
	
	  if answer=='y'
	    loadval=input(' Enter magnitude of load. ');
	    force(i,:)=[nodenum 2 loadval];
	  end
	
	  if answer=='theta'
	    loadval=input(' Enter magnitude of load. ');
	    force(i,:)=[nodenum 3 loadval];
	  end
	
	  answer=input(' Load another node? (y/n) ','s');
	  if answer~='y',break,end
	  plot(node(nodenum,1),node(nodenum,2),'xi')
	  plot(node(nodenum,1),node(nodenum,2),'*b')
	end
	node;
	ncon;
	zero;
	force;
	elseif strcmp(mode,'new')
	disp('new')
	    node=zeros(1000,2);
		ncon=zeros(1000,7);
		zero=zeros(1000,2);
		conm=zeros(1000,2);
		force=zeros(1000,3);
		userdata{1}=node;
		userdata{2}=ncon;
		userdata{3}=conm;
		userdata{4}=zero;
		userdata{5}=force;
		set(gcf,'Userdata',userdata)

	elseif strcmp(mode,'open')
		[filename,pathname] = uigetfile('*_ang*','Select FEA file');

	if filename==0
		break
	end
	
	delete(findobj('tag','opendlg'))
	
	Data=1:64;Data=(Data'*Data)/64;
	msgh=msgbox({'Initiating Load.' '' '' ''},'','custom',Data,hot(64));
	set(msgh,'tag','opendlg')
	mschil=get(msgh,'children');
	delete(findobj(mschil,'style','pushbutton'));
	sh=findobj(get(msgh,'children'),'type','uicontrol');
	%disp('hi there')
	set(sh,'string',['Opening ' filename '.'])
	
	fid=fopen([pathname filename],'rt');
	clear a time phi psi tht htracker7 htracker8 htracker9
	set(sh,'string','Clearing old data')
	
	dfile=dir([pathname filename]);
	set(sh,'string',{['Loading ' filename '.'] 'Please wait.'})
	fprintf('File is %2.1f MB. \n',dfile.bytes/2^20)
	fprintf('Please wait.\n')
	
	
	line = fgetl(fid);
	
	a=fscanf(fid,'%g %g %g %g %g %g %g',[7 inf]);
	a=a';
	if strcmp('T',line(1))~=1
	%	set(sh,'string',{'First line of data file missing'...
	%	'Correcting'})
		a=[str2num(line);a];
		set(sh,'string','Data reconstructed')
	end
	a(:,1);
	time=a(:,1);
	phi=a(:,2);
	psi=a(:,3);
	tht=a(:,4);
	htracker7=a(:,5);
	htracker8=a(:,6);
	htracker9=a(:,7);
	  
	
	
	st=fclose(fid);
	
	set(sh,'string',{[filename ' closed.'] 'Load complete.'})
	%eval(['save ' filename '.mat'])
	%disp([filename '.mat saved.'])
	
	
	if gzipped==1
		set(sh,'string',['Gzipping ' filename ' on your ' computer...
		' system in the background.'])
		if strcmp('MAC2',computer)
			gzip(pathname,filename);
			%figure(msgh)
		elseif isunix
			eval(['!gzip ' pathname filename ' &'])
		end
	end
	
	pause(.2)
	delete(msgh)
	

		
	
	elseif strcmp(mode,'save')
    disp('save')
	file_save_m 
	%path(path,pwd)
	%answer=input(' Save configuration file (Else all will have been in vain)? (y/n) ','s');
	%if answer=='y'
	%  filename=input(' Enter name of configuration file. ','s');
	%  eval(['save ',filename,'.con',' node',' ncon',' zero',' force',' conm']);
	%  answer=input(' Run analysis? (y/n) ','s');
	%  if answer=='y'
	%    vtb8_2(filename);
	%    vtb8_2(node,ncon,zero,force,conm);
	%  end
	%end
	if answer=='y'
	  [filename,pathname]=uiputfile('projectname.con','Save as:');
	  sfilename=size(filename);
	  path(path,pathname)
	  lfilename=sfilename(2);
	  if filename(lfilename-3)=='.'
	    filename=[filename(1:lfilename-4)  '.con'];
	    projectname=filename(1:lfilename-4);
	   else
	    filename=[filename '.con'];
	    projectname=filename(1:lfilename-4);
	  end
	  pathname;
	  sizepath=size(pathname);
	  shortpathname=pathname(1:sizepath(2)-1);
	  lsp=size(findstr(shortpathname,':'));
	  if strcmp(computer,'MAC2') & lsp(1)==0
	    shortpathname=[shortpathname ':'];
	  end
	  %path(shortpathname)
	  cdpath=['cd ' '''' shortpathname '''' ];% Crazy quotes allow spaces
	  %                                       % in directory names.
	  eval(cdpath)
	  eval(['save ',filename,' node',' ncon',' zero',' force',' conm']);
	  answer=input(' Run analysis? (y/n) ','s');
	  if answer=='y'
	    vtb8_2(projectname);
	  end
	end
end
end
