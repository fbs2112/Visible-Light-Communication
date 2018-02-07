%Volterra NLMS Equalyzer using PAM symbols, different evaluation of SNR and of
%the nonlinearity



clear;
clc;
close all;

addpath(['..' filesep 'Equalizer' filesep 'resultsMSE']);
addpath(['.' filesep 'LED Parameters']);

load whiteLED_334-15.mat;

filterFile = 'testDFEVolterraEq';
load([filterFile '.mat']);


%-------------------------Adaptive Filtering Parameters--------------------

M = 4;
numberOfSymbols = 40000;
numberOfBits = log2(M);

blockLength = numberOfSymbols*numberOfBits;
monteCarloLoops = 1000;

volterraFFFlag = 1;
volterraFBFlag = 0;


barGammaVector = 1;
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

if ~volterraFFFlag
    adaptfiltFF = feedforwardLength;
end

if ~volterraFBFlag
    adaptfiltFB = feedbackLength;
end

adapFiltLength = adaptfiltFF + adaptfiltFB;


%-------------------------Adaptive Filtering Parameters--------------------



%-------------------------LED Parameters-----------------------------------
maxLEDVoltage = 3.6; %500 mV
minLEDVoltage = 3;
maxLEDCurrent = 0.03; %500 mA
minLEDCurrent = 0.004; %500 mA

LEDArea = pi*(5e-3)^2;

maxElectricalPower = maxLEDVoltage*maxLEDCurrent;
minElectricalPower = minLEDCurrent*minLEDVoltage;
ISat = ISat;
VB = 2.6; %minimum voltage for current flow 
nLED = n; %LED ideality factor
VT = 0.025; %Thermal voltage


halfAngleLED = deg2rad(15);
luminousIntensityLED = 21375; %milicandela
maxLuminousIntensityLED = 28500;%milicandela

maxCd = 28.5;
minCd = 14.25;

ledLuminousEfficacy = (maxCd - minCd)/(maxElectricalPower - minElectricalPower) ; %this electrical power is evaluated using current and voltage of the linear region of the I-V curve%

fs = 2e6;


Poptical = @(ledLuminousEfficacy,electricalPower,k) (ledLuminousEfficacy.*electricalPower.*LEDArea)./((1 + (ledLuminousEfficacy.*electricalPower.*LEDArea./(maxLuminousIntensityLED.*LEDArea/1000)).^(2*k)).^(1/(2*k)));

%-------------------------LED Parameters-----------------------------------

%-------------------------Photodiode Parameters----------------------------

A = 1e-4; %photodiode area (cm)
d = 10e-2; %distance between LED and photodiode (cm)
R = 0.5;
FOV = deg2rad(25);

%-------------------------Photodiode Parameters----------------------------


%-------------------------Pre Amplifier Parameters-------------------------

transimpedanceGain = 10; %not used

%-------------------------Pre Amplifier Parameters-------------------------


%-------------------------Transmission Parameters--------------------------

kNonLinearity = 2;

theta = 0;
phi = 0;

n = -log(2)/log(cos(halfAngleLED));

H_0 = A/d^2 * (n+1)/(2*pi) * cos(phi)^n * cos(theta) * rectangularPulse(-1,1,theta/FOV);

VDC = 3.25; 
maxAbsoluteValueModulation = 3;

maxModulationIndex = (maxLEDVoltage - VDC)/VDC;

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

    for barGammaIndex = 1:length(barGammaVector)
        equalizerFilter = wFinal(index,:).';

        for SNRIndex = 1:length(SNR)

            for j = 1:monteCarloLoops
                j
                equalizedSignal = zeros(numberOfSymbols,1);

                binaryInputData = randi([0,1],blockLength+100,1);
                binaryInputData = reshape(binaryInputData,[],numberOfBits);
                deciInputData = bi2de(binaryInputData);    
                pilot = real(pammod(deciInputData,2^numberOfBits,0,'gray'));

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

                iLEDOutput = I_V_Fun(filteredVin,VT,nLED,ISat);
                eletricalPowerOutput = filteredVin.*iLEDOutput;

                opticalPowerOutput = Poptical(ledLuminousEfficacy,eletricalPowerOutput,kNonLinearity);

                opticalPowerOutputConvolved = opticalPowerOutput*H_0;

                 n = randn(length(opticalPowerOutputConvolved),1); %noise signal

                receivedCurrentSignal = opticalPowerOutputConvolved*R*A;
                receivedCurrentSignalAC = receivedCurrentSignal - mean(receivedCurrentSignal);
                receivedCurrentSignalPower = receivedCurrentSignalAC'*receivedCurrentSignalAC/length(receivedCurrentSignal);

                powerNoiseAux = n'*n/(length(n));
                powerNoise = (receivedCurrentSignalPower/db2pow(SNR(SNRIndex)));
                n = n.*sqrt(powerNoise/powerNoiseAux);

                receivedVoltageSignalAux = (receivedCurrentSignal + n);
                receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux);
                transimpedanceGain = maxAbsoluteValueModulation/max(receivedVoltageSignalAux);
                receivedVoltageSignal =  receivedVoltageSignalAux*sqrt(var(pilot)/var(receivedVoltageSignalAux));

                unbiasedReceivedVoltageSignal = receivedVoltageSignal - VDC;

                xAux = [zeros(N-1,1);receivedVoltageSignal];

                outputFB = zeros(length(pilot),1);
                outputFF = zeros(length(pilot),1);
                inputFB = zeros(feedbackLength,length(pilot));
                
                for k = feedforwardLength:length(pilot)

                    xFlip = xAux(k:-1:k-feedforwardLength+1);

                     if volterraFFFlag

                        aux = zeros((feedforwardLength^2+feedforwardLength)/2,1);

                        for lIndex = 1:length(l1FF)
                            aux(lIndex,1) = xFlip(l1FF(lIndex),1)*(xFlip(l2FF(lIndex),1));
                        end
                        xConc = [xFlip(:,1);aux];
                    else
                        xConc = xFlip(:,1);
                     end


                     if ~volterraFFFlag && ~volterraFBFlag 
                        xConc = x(:,k);
                     end
                     
                    outputFF(k) =  (equalizerFilter(1:adaptfiltFF))'*xConc;
                    equalizedSignal(k,1) = pamHardThreshold(outputFF(k) + outputFB(k-1)); %MEXER AQUI
                   
                    inputFB(:,k) = equalizedSignal(k:-1:k-feedbackLength + 1); 


                    if volterraFBFlag
                        aux = zeros((feedbackLength^2+feedbackLength)/2,1);
                        for lIndex = 1:length(l1FB)
                            aux(lIndex,1) = inputFB(l1FB(lIndex),k)*(inputFB(l2FB(lIndex),k));
                        end

                        yHatConc = [inputFB(:,k);aux];
                    else
                        yHatConc = inputFB(:,k);
                    end

                    if ~volterraFFFlag && ~volterraFBFlag 
                        yHatConc = inputFB(:,k);
                    end

                    outputFB(k) = (equalizerFilter(adaptfiltFF+1:end))'*yHatConc;

                end
                
               [corr,lags] = xcorr(equalizedSignal,xAux(N:end));
               [~,idx] = max(abs(corr));
               delay = abs(lags(idx));

               decDemodSignal = pamdemod(equalizedSignal,2^numberOfBits,0,'gray');
               binaryOutputData = de2bi(decDemodSignal,numberOfBits);
               berAux(j) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;               
            end
            
            
            ber(index,SNRIndex) = mean(berAux);
           
        end
   end

end

save(['.' filesep 'resultsBER' filesep filterFile '_BER.mat'],'ber','SNR');
rmpath(['..' filesep 'Equalizer' filesep 'resultsMSE']);
rmpath(['.' filesep 'LED Parameters']);
