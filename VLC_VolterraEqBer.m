%%
%Volterra NLMS Equalyzer using PAM symbols, different evaluation of SNR and of
%the nonlinearity



clear;
clc;
close all;

addpath(['..' filesep 'Equalyzer' filesep 'resultsMSE_VLC']);
addpath(['.' filesep 'LED Parameters']);

load whiteLED_334-15.mat;

load results30.mat;

% load channel01.mat;



% R = 0.56;
% R = 1;
% 
% 
% b0 = 1;
% b1 = 0.5;
% b2 = 0.05;

%-------------------------Adaptive Filtering Parameters--------------------
M = 4;
numberOfSymbols = 40000;
numberOfBits = log2(M);

blockLength = numberOfSymbols*numberOfBits;
monteCarloLoops = 1000;

% auxMatrix = triu(ones(N));
% [l1,l2] = find(auxMatrix);
% 
% delayVector = N/2;


% noisePower = 100;
% % barGamma = 4*sqrt(5*noisePower);
% 
% barGamma = 0;
% 
barGammaVector = 1;
N = 12;


%-------------------------Adaptive Filtering Parameters--------------------



adapFiltLength = (N^2+N)/2 + N;

% adapFiltLength = N;


auxMatrix = triu(ones(N));
[l1,l2] = find(auxMatrix);




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

transimpedanceGain = 10;

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

SNR = 0:5:30;

berAux = zeros(monteCarloLoops,1);
ber = zeros(length(SNR),length(modulationIndexVector));




for index = 1:length(modulationIndexVector)
    
    modulationIndex = modulationIndexVector(index);
    
    if modulationIndex > maxModulationIndex
        warning('Modulation Index may cause undesired nonlinear effects')
    end

    maxVoltage = VDC*(1+modulationIndex);
    deltaV = maxVoltage - VDC;


    VoltageConstant = modulationIndex*maxVoltage/((1+modulationIndex)*maxAbsoluteValueModulation);


    for barGammaIndex = 1:length(barGammaVector)
        equalyzerFilter = squeeze(wFinal(index,1,1,:));


        for SNRIndex = 1:length(SNR)


            for j = 1:monteCarloLoops
                j
                equalyzedSignal = zeros(numberOfSymbols,1);


                binaryInputData = randi([0,1],blockLength+100,1);
                binaryInputData = reshape(binaryInputData,[],numberOfBits);
%                 deciInputData = randi([0,M-1],2000,1);
                deciInputData = bi2de(binaryInputData);    
                pilot = real(pammod(deciInputData,2^numberOfBits,0,'gray'));

%                     input = randi([0,2^numberOfBits-1],maxRuns*2,1);
%                     pilot = pammod(input,2^numberOfBits,0,'gray');

%                 Vin = pilot*VoltageConstant + VDC; %Using symbols to modulate voltage
                Vin = pilot;
                convLength = length(Vin) + 1000 -1;
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



%                 opticalPowerOutput = eletrical2OpticalGain*eletricalPowerOutput;

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
                powerNoise = (receivedCurrentSignalPower/db2pow(SNR(SNRIndex)));
                n = n.*sqrt(powerNoise/powerNoiseAux);
    %                 
                receivedVoltageSignalAux = (receivedCurrentSignal + n);
                receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux);
                transimpedanceGain = maxAbsoluteValueModulation/max(receivedVoltageSignalAux);
                receivedVoltageSignal =  receivedVoltageSignalAux*sqrt(var(pilot)/var(receivedVoltageSignalAux));



                unbiasedReceivedVoltageSignal = receivedVoltageSignal - VDC;

                xAux = [zeros(N-1,1);receivedVoltageSignal];



                for k = N:length(pilot)

                    xFlip = xAux(k:-1:k-N+1);
                    
                    
%                     xAP = zeros(N,L+1);
% 
%                     for l = 0:L
%                         xAP(:,l+1) = xAux(k-l:-1:k-N+1-l);
%                     end
% 
% %                     xAP = flipud(xAux);

%                     xTDLConc = zeros(adapFiltLength,1);

                   

                    xTDLAux = zeros((N*N+N)/2,1);

                    for lIndex = 1:length(l1)
                        xTDLAux(lIndex,1) = xFlip(l1(lIndex),1)*(xFlip(l2(lIndex),1));

%                                 xTDLAux(counterAux,1) = xTDL(l1,k+l3-1)*xTDL(l2,k+l3-1);
                    end


                    xTDLConc = [xFlip(:,1);xTDLAux];
        %                 xTDLConc(:,l3+k-L-1) = [xTDL(:,k+l3-1);xTDLAux];
                    



                    equalyzedSignal(k,1) =  (equalyzerFilter)'*xTDLConc;


                end
%                 equalyzedSignal2 = (equalyzedSignal - VDC)/VoltageConstant;
                
               [corr,lags] = xcorr(equalyzedSignal,xAux(N:end));
               [~,idx] = max(abs(corr));
               delay = abs(lags(idx));

            %     equalyzedSignal = equalyzedSignal(1:1000,:);
%                equalyzedSignal = 
               decDemodSignal = pamdemod(equalyzedSignal,2^numberOfBits,0,'gray');
               
%                berAux(j) = sum(sum(abs(equalyzedSignal(delay+1:end,:) - Vin(1:end-delay,:))))./length(equalyzedSignal);
               binaryOutputData = de2bi(decDemodSignal,numberOfBits);


               berAux(j) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
               
            end
            
            
            ber(SNRIndex,index) = mean(berAux);
           

        end
   

    end

end



save(['.' filesep 'resultsBER' filesep 'resultsBER22.mat'],'ber','SNR');



rmpath(['..' filesep 'Equalyzer' filesep 'resultsMSE_VLC']);
rmpath(['.' filesep 'LED Parameters']);




