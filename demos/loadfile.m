%LOADFILE loads a file using the GUI open dialog box.
% The selected MAT-file will be loaded.
%
% See also ADDPATH, DELETEFILE, EDITFILE, FILEBAR

% Copyright Joseph C. Slater 1994

[filename,pathname]=uigetfile('*','Select M-File to Load');
destr=['load ' '''' [pathname filename] '''' ' -mat'];
if exist([pathname filename])==1 | exist([pathname filename])==2
  eval(destr);
  disp(['File ' filename ' loaded.'])
 else
  disp(['File ' filename ' not found. No action taken'])
end
clear pathname filename
