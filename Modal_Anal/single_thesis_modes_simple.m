% thesis
clear all
close all

%C:\WINDOWS\Desktop\SINGLE\L5_short_tubes.trf\L5_S02_ASCIIdata\
Hxy = uigetfile('*.txt','Pick a Hxy (transfer function) file');

% get transfer function from file
% ASCII file setup is {Hz Real Imag}
% Channel 2 - x dir
% Channel 3 - y dir
% Channel 4 - z dir
% Channel 5 - velocity from vibrometer
% Channel 8 - base accelerometer


%data =  dlmread(Hxy,'\t',8,1);
load dat
f = data(:,1);  % in Hz

x = data(:,6);
y = i*data(:,7); % z dir in g/g


TF = x+y;

[z,nf,u]=mmcf(f,TF)
%%% need mmcf.p to run this script %%%

%[z,nf,u]=mmcf(f,TF,Fmin,Fmax,nmodes) Curve fit to MDOF FRF.
% f is the frequency vector in Hz. 
%
% TF is the complex transfer functions, each FRF in a column.
% z and nf are the damping ratio and natural frequency (Hz)
% u is the mode shape.
%
% If Fmin and Fmax are not entered, the min and max values in
% f are used.
%
% If nmodes is not included, it is estimated automatically.

