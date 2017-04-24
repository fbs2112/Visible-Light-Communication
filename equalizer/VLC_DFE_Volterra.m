%%
%Volterra NLMS Equalyzer using PAM symbols, different evaluation of SNR and of
%the nonlinearity



clear;
clc;
close all;

addpath(['..' filesep '..' filesep 'Channel results']);
addpath(['..' filesep 'VLC Simulator' filesep]);
addpath(['..' filesep 'VLC Simulator' filesep 'LED Parameters']);

load whiteLED_334-15.mat;

% load channel01.mat;


% R = 0.56;
% R = 1;
% 
% 
% b0 = 1;
% b1 = 0.5;
% b2 = 0.05;

%-------------------------Adaptive Filtering Parameters--------------------
numberOfBits = 2;
maxRuns = 15000;
maxIt = 1000;
gamma = 1e-12;
SNR = 30;
mu = 0.8;

volterraFFFlag = 1;
volterraFBFlag = 0;

feedforwardLength = 12;
feedbackLength = 12;

N = feedforwardLength;

adaptfiltFF = (feedforwardLength^2+feedforwardLength)/2 + feedforwardLength;
adaptfiltFB = (feedbackLength^2+feedbackLength)/2 + feedbackLength;

adaptfilt = adaptfiltFF + adaptfiltFB;

auxMatrix = triu(ones(feedforwardLength));
[l1FF,l2FF] = find(auxMatrix);

auxMatrix = triu(ones(feedbackLength));
[l1FB,l2FB] = find(auxMatrix);
% delayVector = 1:15;
delayVector = round((feedforwardLength)/2);

delayVector = 14;
% delayVector = 10;


if ~volterraFFFlag
    adaptfiltFF = feedforwardLength;
end

if ~volterraFBFlag
    adaptfiltFB = feedbackLength;
end

    
adapFiltLength = adaptfiltFF + adaptfiltFB;


noisePower = 100;
% barGamma = 4*sqrt(5*noisePower);

barGamma = 0;

barGammaVector = 1;

%-------------------------Adaptive Filtering Parameters--------------------



%-------------------------LED Parameters-----------------------------------
maxLEDVoltage = 3.6; %500 mV
minLEDVoltage = 3;
maxLEDCurrent = 0.03; %500 mA
minLEDCurrent = 0.004; %500 mA

maxElectricalPower = maxLEDVoltage*maxLEDCurrent;
minElectricalPower = minLEDCurrent*minLEDVoltage;
% TOV = 0.2; 
% eletrical2OpticalGain = 1; %eletrical to optical gain imposed by the LED

ISat = ISat;
VB = 2.6; %minimum voltage for current flow 
nLED = n; %LED ideality factor
VT = 0.025; %Thermal voltage


halfAngleLED = deg2rad(15);
luminousIntensityLED = 21375; %milicandela
maxLuminousIntensityLED = 28500;%milicandela

% opticalPower = luminousIntensityLED*2*pi*(1-cos(halfAngleLED))/1000;

% ledLuminousEfficacy = opticalPower/(3.2*10e-3); %this electrical power is evaluated using current and voltage of the linear region of the I-V curve%
maxCd = 28.5;
minCd = 14.25;




% ledLuminousEfficacy = opticalPower/(3.2*10e-3); %this electrical power is evaluated using current and voltage of the linear region of the I-V curve%
ledLuminousEfficacy = (maxCd - minCd)/(maxElectricalPower - minElectricalPower) ; %this electrical power is evaluated using current and voltage of the linear region of the I-V curve%


fs = 2e6;

% f = fs/2*linspace(0,1,1000) *2*pi;
% 
% w = [-fliplr(f(2:end-1)) f];
% 
% LEDResp = freqRespLED(w);




%  NFFT = 2^nextpow2(length(LEDResp)); % Next power of 2 from length of y
% %                 Y = fft(noise,length(LEDResp))/length(LEDResp);
%                 f = fs/2*linspace(0,1,length(LEDResp));
% 
%                 % Plot single-sided amplitude spectrum.
%                 figure;
%                 plot(f,20*log10((LEDResp(1:length(LEDResp)))) )
% 


Poptical = @(ledLuminousEfficacy,electricalPower,k) (ledLuminousEfficacy.*electricalPower)./((1 + (ledLuminousEfficacy.*electricalPower./(maxLuminousIntensityLED/1000)).^(2*k)).^(1/(2*k)));

%-------------------------LED Parameters-----------------------------------




%-------------------------Photodiode Parameters----------------------------

A = 1e-4; %photodiode area (cm)
d = 10e-2; %distance between LED and photodiode (cm)
R = 0.5;

FOV = deg2rad(25);

%-------------------------Photodiode Parameters----------------------------


%-------------------------Pre Amplifier Parameters-------------------------

transimpedanceGain = 1;

%-------------------------Pre Amplifier Parameters-------------------------


%-------------------------Transmission Parameters--------------------------

kNonLinearity = 2;


% bitRate = 1e6; %1 Mb/s

theta = 0;
phi = 0;

n = -log(2)/log(cos(halfAngleLED));

H_0 = A/d^2 * (n+1)/(2*pi) * cos(phi)^n * cos(theta) * rectangularPulse(-1,1,theta/FOV);


VDC = 3.25;
maxAbsoluteValueModulation = 3;

maxModulationIndex = (maxLEDVoltage - VDC)/VDC;
% modulationIndexVector = 0.01:0.02:maxModulationIndex;
modulationIndexVector = [0.05 0.075 0.1];


%-------------------------Transmission Parameters--------------------------




for index = 1:length(modulationIndexVector)
    modulationIndex = modulationIndexVector(index);
    
    if modulationIndex > maxModulationIndex
        warning('Modulation Index may cause undesired nonlinear effects')
    end

    maxVoltage = VDC*(1+modulationIndex);
    deltaV = maxVoltage - VDC;


    

    for barGammaIndex = 1:length(barGammaVector)
        count = zeros(maxIt,1);

    %     count = zeros(maxIt,length(barGammaVector));
        for delay = 1:length(delayVector)

            for L = 0:0

                u = zeros(L+1,1);
                u(1) = 1;


                w2 = zeros(adapFiltLength,maxRuns,maxIt);
                for j = 1:maxIt
                    j


                    input = randi([0,2^numberOfBits-1],maxRuns*2,1);
%                     blockLength = 2*maxRuns;
%                     binaryInputData = randi([0,1],blockLength+100,1);
%                     binaryInputData = reshape(binaryInputData,[],numberOfBits);
% %                   deciInputData = randi([0,M-1],2000,1);
%                     input = bi2de(binaryInputData);    
                    pilot = real(pammod(input,2^numberOfBits,0,'gray'));

%                     Vin = pilot*VoltageConstant + VDC; %Using symbols to modulate voltage
                    Vin = pilot;

                    convLength = maxRuns*2 + 1000 -1;
                    NFFT = 2^nextpow2(convLength);

                    VinFreq = fft(Vin,NFFT);

                    f = fs/2*linspace(0,1,NFFT/2 + 1)  *2*pi;

                    w = [-fliplr(f(2:end-1)) f];

                    LEDResp = freqRespLED(w);

                    filteredVinAux = real(ifft(VinFreq.*fftshift(LEDResp))); 

                    filteredVin = filteredVinAux(1:length(Vin));
                    
                    VoltageConstant = modulationIndex*maxVoltage/((1+modulationIndex)*max(filteredVin));
                
                    Vin = filteredVin*VoltageConstant + VDC;
                    filteredVin = Vin;



    %                 iLEDOutput = ledModel(I_V_Fun,Vin,maxLEDVoltage,kNonLinearity);

                    iLEDOutput = I_V_Fun(filteredVin,VT,nLED,ISat);                
    %                 iLEDOutput = Vin;

    %                 iLEDOutput = 1;

                    eletricalPowerOutput = filteredVin.*iLEDOutput;



    %                 opticalPowerOutput = ledLuminousEfficacy*eletricalPowerOutput;

                    opticalPowerOutput = Poptical(ledLuminousEfficacy,eletricalPowerOutput,kNonLinearity);

                    opticalPowerOutputConvolved = opticalPowerOutput*H_0;






    %                 NFFT = 2^nextpow2(length(LEDResp)); % Next power of 2 from length of y
    % %                 Y = fft(noise,length(LEDResp))/length(LEDResp);
    %                 f = fs/2*linspace(0,1,length(LEDResp));
    % 
    %                 % Plot single-sided amplitude spectrum.
    %                 figure;
    %                 plot(f,20*log10(ifftshift(LEDResp(1:length(LEDResp)))) )
    %                 
    %                 figure
    %                 
    %                 
    %                 plot(f,20*log10(abs(VinFreq)))
    %                 
    %                 figure
    %                 aux = VinFreq.*ifftshift(LEDResp);
    %                 
    %                 plot(f,20*log10(abs(aux)))
    %                 



    %                 Pout = mean(opticalPowerOutput);


    %                 varNoise = (R^2*H_0^2*Pout^2)/(bitRate*db2pow(SNR));


    %                 n = n.*sqrt(varNoise/powerNoiseAux);
                    n = randn(length(opticalPowerOutputConvolved),1); %noise signal

                    receivedCurrentSignal = opticalPowerOutputConvolved*R*A;
                    receivedCurrentSignalAC = receivedCurrentSignal - mean(receivedCurrentSignal);
                    receivedCurrentSignalPower = receivedCurrentSignalAC'*receivedCurrentSignalAC/length(receivedCurrentSignal);

                    powerNoiseAux = n'*n/(length(n));
                    powerNoise = (receivedCurrentSignalPower/db2pow(SNR));
                    n = n.*sqrt(powerNoise/powerNoiseAux);

    %                 

                    receivedVoltageSignalAux = (receivedCurrentSignal + n);
                    receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux);
                    transimpedanceGain = maxAbsoluteValueModulation/max(receivedVoltageSignalAux);
                    receivedVoltageSignal =  receivedVoltageSignalAux*sqrt(var(pilot)/var(receivedVoltageSignalAux));
                    
                    unbiasedReceivedVoltageSignal = receivedVoltageSignal - VDC;

                    xAux = [zeros(N-1,1);receivedVoltageSignal];

                    w = zeros(adapFiltLength,maxRuns);

        %                 w = zeros(N,maxRuns); 
%                     pilot = Vin;


                    for k = (adapFiltLength + delayVector(delay) + L + 10):maxRuns

                        x(:,k) = xAux(k:-1:k-feedforwardLength+1);

                        yHat(:,k) = (pilot(-delayVector(delay) + k + 1 -1:-1:-delayVector(delay) + k + 1 - feedbackLength - 1 + 1));

                        if volterraFFFlag

                            aux = zeros((feedforwardLength^2+feedforwardLength)/2,1);

                            for lIndex = 1:length(l1FF)
                                aux(lIndex,1) = x(l1FF(lIndex),k)*(x(l2FF(lIndex),k));
                            end
                            xConc = [x(:,k);aux];
                        else
                            xConc = x(:,k);
                        end


                        if volterraFBFlag
                            aux = zeros((feedbackLength^2+feedbackLength)/2,1);
                            for lIndex = 1:length(l1FB)
                                aux(lIndex,1) = yHat(l1FB(lIndex),k)*(yHat(l2FB(lIndex),k));
                            end

                            yHatConc = [yHat(:,k);aux];
                        else
                            yHatConc = yHat(:,k);
                        end

                        if ~volterraFFFlag && ~volterraFBFlag 
                            xConc = x(:,k);
                            yHatConc = yHat(:,k);
                        end

                        z = [xConc;yHatConc];

    %                     xAP = zeros(N,L+1);
    % 
    %                     for l = 0:L
    %                         xAP(:,l+1) = xAux(k-l:-1:k-N+1-l);
    %                     end
    % 
    % %                     xAP = flipud(xAux);
    % 
    %                     xTDLConc = zeros(adapFiltLength,L+1);
    % 
    %                     for l3 = 1:L+1
    %                         xTDLAux = zeros((N*N+N)/2,1);
    % 
    % 
    % 
    %                         for lIndex = 1:length(l1)
    %                             xTDLAux(lIndex,1) = xAP(l1(lIndex),l3)*(xAP(l2(lIndex),l3));
    % 
    %     %                                 xTDLAux(counterAux,1) = xTDL(l1,k+l3-1)*xTDL(l2,k+l3-1);
    %                         end
    % 
    % 
    %                         xTDLConc(:,l3) = [xAP(:,l3);xTDLAux];
    %         %                 xTDLConc(:,l3+k-L-1) = [xTDL(:,k+l3-1);xTDLAux];
    %                     end

                        if k == 15000
                           k; 
                        end

                        d(k) = (pilot(-delayVector(delay) + k + 1)); 


                        e(k) = d(k) - w(:,k)'*z;

                        w(:,k+1) = w(:,k) + mu*z*((z'*z+gamma*eye(L+1))\eye(L+1))*conj(e(k))*u;

    %                   
                    end
                    w2(:,:,j) = (w(:,1:maxRuns));
                    e2(:,j) = abs(e).^2;
                end

                meanCount(barGammaIndex) = mean(count);

    %             count = zeros(maxIt,1);

                w3 = mean(w2,3);
                wFinal(index,barGammaIndex,delay,:,L+1) = w3(:,end);
                
%                 wFinal = w3;


                e3(index,barGammaIndex,delay,:,L+1) = mean(e2,2);

            end
            % save(['.' filesep 'results' filesep 'results07.mat'],'wFinal','e3','meanCount');




        %     end

        end

    end

end




% for i = 1:length(delayVector)
%     figure
% plot(10*log10((e3(i,:))))
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% end


% for i = 1:adapFiltLength+10
%     plot(10*log10((e3(i,:,1))))
%     xlabel('Iterations','interpreter','latex');
%     ylabel('MSE (dB)','interpreter','latex');
%     hold on;
% end


save(['.' filesep 'resultsMSE_VLC' filesep 'results32.mat'],'wFinal','e3','meanCount');



rmpath(['..' filesep 'VLC Simulator' filesep]);
rmpath(['..' filesep '..' filesep 'Channel results']);
rmpath(['..' filesep 'VLC Simulator' filesep 'LED Parameters']);




