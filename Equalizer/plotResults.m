clear;
clc;
close all;


addpath(['.' filesep 'resultsMSE']);

load testLinEq.mat;

for i = 1:size(e3,1)
    x = find(squeeze(e3(i,:)),1);
    plot(1:length(e3)-x,10*log10(squeeze(e3(i,x:end-1))));

    hold on;

end
xlabel('Iterations','interpreter','latex');
ylabel('MSE (dB)','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
xlim([0 500]);
ylim([-20 4]);


figure
load testVolterraEq.mat;


for i = 1:size(e3,1)
    x = find(squeeze(e3(i,:)),1);
    plot(1:length(e3)-x,10*log10(squeeze(e3(i,x:end-1))));

    hold on;

end
xlabel('Iterations','interpreter','latex');
ylabel('MSE (dB)','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
% xlim([0 500]);
ylim([-20 4]);


figure
load testDFEEq.mat;


for i = 1:size(e3,1)
    x = find(squeeze(e3(i,:)),1);
    plot(1:length(e3)-x,10*log10(squeeze(e3(i,x:end-1))));

    hold on;

end
xlabel('Iterations','interpreter','latex');
ylabel('MSE (dB)','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
xlim([0 500]);
ylim([-20 4]);



figure
load testDFEVolterraEq.mat;


for i = 1:size(e3,1)
    x = find(squeeze(e3(i,:)),1);
    plot(1:length(e3)-x,10*log10(squeeze(e3(i,x:end-1))));

    hold on;

end
xlabel('Iterations','interpreter','latex');
ylabel('MSE (dB)','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
% xlim([0 500]);
ylim([-20 4]);


