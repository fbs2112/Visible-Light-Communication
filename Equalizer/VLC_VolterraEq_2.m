%%
%Volterra NLMS Equalyzer using PAM symbols, different evaluation of SNR and of
%the nonlinearity

clear;
clc;
close all;

addpath(['..' filesep 'VLC_Simulator' filesep]);
addpath(['..' filesep 'VLC_Simulator' filesep 'LED Parameters']);

load whiteLED_334-15.mat;

%-------------------------Adaptive Filtering Parameters--------------------
numberOfBits = 2;
N = 12;
maxRuns = 15000;
maxIt = 100;
mu = 0.8;
gamma = 1e-12;
SNR = 30;


adapFiltLength = (N^2+N)/2 + N;

auxMatrix = triu(ones(N));
[l1,l2] = find(auxMatrix);


delayinSamples = 13;

%-------------------------Adaptive Filtering Parameters--------------------


%-------------------------LED Parameters-----------------------------------

Poptical = @(ledLuminousEfficacy,electricalPower,k) (ledLuminousEfficacy.*electricalPower)./((1 + (ledLuminousEfficacy.*electricalPower./(maxLuminousIntensityLED/1000)).^(2*k)).^(1/(2*k)));

%-------------------------LED Parameters-----------------------------------


%-------------------------Photodiode Parameters----------------------------

A = 1e-4; %photodiode area (cm)
d = 10e-2; %distance between LED and photodiode (cm)
R = 0.5;
FOV = deg2rad(25);

%-------------------------Photodiode Parameters----------------------------


%-------------------------Pre Amplifier Parameters-------------------------

%-------------------------Pre Amplifier Parameters-------------------------


%-------------------------Transmission Parameters--------------------------

kNonLinearity = 2;

LEDfreqRespPoints = 1000;

fs = 2e6;

theta = 0;
phi = 0;

H_0 = A/d^2 * (n+1)/(2*pi) * cos(phi)^n * cos(theta) * rectangularPulse(-1,1,theta/FOV);

VDC = 3.25;
maxAbsoluteValueModulation = 3;

maxModulationIndex = (maxLEDVoltage - VDC)/VDC;
modulationIndexVector = [0.05 0.075 0.1];


%-------------------------Transmission Parameters--------------------------


wFinal = zeros(length(modulationIndexVector),adapFiltLength);
e3 = zeros(length(modulationIndexVector),maxRuns + delayinSamples + 1);

for index = 1:length(modulationIndexVector)
    modulationIndex = modulationIndexVector(index);
    
    if modulationIndex > maxModulationIndex
        warning('Modulation Index may cause undesired nonlinear effects')
    end

    maxVoltage = VDC*(1+modulationIndex);
    deltaV = maxVoltage - VDC;

   
    w2 = zeros(adapFiltLength,maxRuns,maxIt);
    e2 = zeros(maxRuns + delayinSamples + 1,maxIt);
    
    for j = 1:maxIt

        input = randi([0,2^numberOfBits-1],maxRuns*2,1);
        pilot = real(pammod(input,2^numberOfBits,0,'gray'));

        convLength = length(pilot) + LEDfreqRespPoints -1;
        NFFT = 2^nextpow2(convLength);

        pilotFreq = fft(pilot,NFFT);

        f = fs/2*linspace(0,1,NFFT/2 + 1)*2*pi;

        w = [-fliplr(f(2:end-1)) f];

        LEDResp = freqRespLED(w);

        filteredVinAux = real(ifft(pilotFreq.*fftshift(LEDResp))); 

        filteredVin = filteredVinAux(1:length(pilot));

        VoltageConstant = modulationIndex*maxVoltage/((1+modulationIndex)*max(filteredVin));

        filteredVin = filteredVin*VoltageConstant + VDC;
            
        iLEDOutput = I_V_Fun(filteredVin,VT,nLED,ISat);

        eletricalPowerOutput = filteredVin.*iLEDOutput;

        opticalPowerOutput = Poptical(ledLuminousEfficacy,eletricalPowerOutput,kNonLinearity);

        opticalPowerOutputConvolved = opticalPowerOutput*H_0;

        n = randn(length(opticalPowerOutputConvolved),1); %noise signal

        receivedCurrentSignal = opticalPowerOutputConvolved*R*A;
        receivedCurrentSignalAC = receivedCurrentSignal - mean(receivedCurrentSignal);
        receivedCurrentSignalPower = receivedCurrentSignalAC'*receivedCurrentSignalAC/length(receivedCurrentSignal);

        powerNoiseAux = n'*n/(length(n));
        powerNoise = (receivedCurrentSignalPower/db2pow(SNR));
        n = n.*sqrt(powerNoise/powerNoiseAux);

        receivedVoltageSignalAux = (receivedCurrentSignal + n);
        receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux);
        receivedVoltageSignal =  receivedVoltageSignalAux*sqrt(var(pilot)/var(receivedVoltageSignalAux));


        xAux = [zeros(N-1,1);receivedVoltageSignal];

        w = zeros(adapFiltLength,maxRuns);

        d = zeros(maxRuns + delayinSamples + 1,1);
        e = zeros(maxRuns + delayinSamples + 1,1);
        
        for k = delayinSamples + 1:maxRuns + delayinSamples + 1

            x = xAux(k:-1:k-N+1);

            xTDLAux = zeros((N*N+N)/2,1);


            for lIndex = 1:length(l1)
                xTDLAux(lIndex,1) = x(l1(lIndex),1)*(x(l2(lIndex),1));
            end
            
            
            xConc = [x;xTDLAux];


            d(k) = (pilot(-delayinSamples + k + 1)); 


            e(k) = d(k) - w(:,k)'*xConc;

            w(:,k+1) = w(:,k) + mu*xConc*((xConc'*xConc+gamma*eye(1))\eye(1))*conj(e(k));  


        end
        w2(:,:,j) = conj(w(:,1:maxRuns));
        e2(:,j) = abs(e).^2;
    end

    w3 = mean(w2,3);
    wFinal(index,:) = w3(:,end);

    e3(index,:) = mean(e2,2);

end


save(['.' filesep 'resultsMSE' filesep 'testVolterraEq.mat'],'wFinal','e3');

rmpath(['..' filesep 'VLC_Simulator' filesep]);
rmpath(['..' filesep 'VLC_Simulator' filesep 'LED Parameters']);
