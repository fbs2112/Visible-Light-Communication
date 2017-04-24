%This script evaluates through an least-square fitting the saturation 
%current and ideality factor of the LED using points of its I-V curve

clear;
close all;
clc;


VT = 25e-3; %thermal voltage


IF = [2 4 6 8 10 12 14 16 18 20] *1e-3; %current

VF = [2.9 3 3.09 3.15 3.2 3.25 3.3 3.39 3.45 3.5]; %voltage

b = log(IF);

m = [ones(length(IF),1) ((VF)/VT).'];

x = (m.'*m)\eye(2)*m.'*b.';


ISat = exp(x(1));

n = 1/x(2);

v = 0:0.01:3.5;

aux = v > 0;

 I = @(ISat,n,VF) (ISat.*exp((VF)./(n.*VT))  - ISat).*aux;
 
 save(['whiteLED_334-15_Param.mat'],'ISat','n');
 


plot(v,I(ISat,n,v))