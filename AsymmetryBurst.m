%--------------------------------------------------------------------------
% Authors: Brandon S Coventry, Claudia M Krogmeier
% Date: 08/03/21
% Purpose: Implement the Asymmetry Burst Analyses in Allen and Cohen 2010
% Inputs: eegsig1 - First electrode for asymmetry calculations
%         eegsig2 - Second electrode for asymmetry calculations
%         time - Time vector in seconds. Should map samples in eegsig1/2 to
%         time.
%         fs - sample rate of EEG acquisition
%         Order matters for eegsig1/2. Asymmetry calculations are performed
%         as follows: log_e(eegsig1)-log_e(eegsig2). For this paper,
%         eegsig1 = F4, eegsig2 = F3
% Outputs:
% Revision Hist: 03/08/2022
%                - Corrected an error in code, streamlined processing
%                pipeline. 
%--------------------------------------------------------------------------
function [asymmetry,waveDecomp,totPeakCount,highPeakCount,lowPeakCount,powersHigh,powersLow,mPowerHigh,mPowerLow,phasesHigh,phasesLow] = AsymmetryBurst(eegsig1,eegsig2,time,fs)
% Begin by setting constants in algorithm
thresholdPercentHigh = 90; %Percentile of bursts to keep. 
thresholdPercentLow = 100-thresholdPercentHigh;   %Lower threshold to keep
burstWin = 0.01;                   %Nominally 10ms
PowerCalculationWindow = 0.05;     %Window over which to calculate average powers.
% Calculate Asymmetry bursts
analyticSig1 = hilbert(eegsig1);          %Calculate the analytic signal with a Hilbert transform. Result is a complex number
analyticSig2 = hilbert(eegsig2);
InstantPower1 = log(abs(analyticSig1).^2);     %Calculates the instantaneous power of the time series from the analytic signal.
InstantPower2 = log(abs(analyticSig2).^2);
asymmetry = InstantPower1-InstantPower2;
% Calculate bursts using Percentiles
lowKeep = prctile(real(asymmetry),thresholdPercentLow);       %Get the Nth percentile (based on thresholdPercentLow/High) 
highKeep = prctile(real(asymmetry),thresholdPercentHigh);
lowIDX = find(real(asymmetry)<lowKeep);
highIDX = find(real(asymmetry)>highKeep);
% Account for sustained burst activity 
[highPeaks,highPeakLocs] = findpeaks(asymmetry,fs,'MinPeakDistance',burstWin,'MinPeakHeight',highKeep);
[lowPeaks,lowPeakLocs] = findpeaks(-1*asymmetry,fs,'MinPeakDistance',burstWin,'MinPeakHeight',-1*lowKeep);
highPeakCount = length(highPeaks);
lowPeakCount = length(lowPeaks);
totPeakCount = highPeakCount+lowPeakCount;
%Find refined indicies
highIDXList = zeros(1,length(highPeakLocs));
lowIDXList = zeros(1,length(lowPeakLocs));
for ck = 1:length(highPeakLocs)
    highIDXList(ck) = find(time==highPeakLocs(ck));
end
for bc = 1:length(lowPeakLocs)
    lowIDXList(bc) = find(time==lowPeakLocs(bc));
end
% Do wavelet decomposition
waveFreqs = logspace(0.8975,1.115);       %Alpha band 8-13 Hz.
waveDecomp = zeros(length(waveFreqs),length(asymmetry));
for bc = 1:length(waveFreqs)
    wavlet = calculateMorlet(1:512,waveFreqs(bc));
    waveDecomp(bc,:) = conv(asymmetry,wavlet,'same');
end
%Do power and phase calculations. 
powersHigh = cell(1,length(highPeakLocs));
mPowerHigh = cell(1,length(highPeakLocs));
phasesHigh = cell(1,length(highPeakLocs));
powersLow = cell(1,length(lowPeakLocs));
mPowerLow = cell(1,length(lowPeakLocs));
phasesLow = cell(1,length(lowPeakLocs));
%Calculate from high powers first
for ck = 1:length(highIDXList)
    curLoc = highIDXList(ck);
    if curLoc > round(fs*PowerCalculationWindow/2)
        if curLoc<length(asymmetry)-round(fs*PowerCalculationWindow/2)
            win = curLoc-round(fs*PowerCalculationWindow/2):curLoc+round(fs*PowerCalculationWindow/2)-1;
        else
            win = curLoc-round(fs*PowerCalculationWindow/2)-(length(asymmetry)-curLoc-1):length(asymmetry);
        end
    else
        win = 1:curLoc+(round(fs*PowerCalculationWindow/2))+(round(fs*PowerCalculationWindow/2)-curLoc);
    end
%     if curLoc > round(fs*PowerCalculationWindow/2)
%         if curLoc<length(asymmetry)-round(fs*PowerCalculationWindow/2)
%             win = curLoc-round(fs*PowerCalculationWindow/2):curLoc+round(fs*PowerCalculationWindow/2)-1;
%         else
%             win = curLoc-round(fs*PowerCalculationWindow/2)-(length(asymmetry)-curLoc-1):length(asymmetry);
%         end
%     else
%         win = 1:curLoc+(round(fs*PowerCalculationWindow/2)-1)+(curLoc-round(fs*PowerCalculationWindow/2)-1);
%     end
    pow = zeros(length(waveFreqs),length(win));
    mpow = zeros(1,length(waveFreqs));
    phaseList = zeros(length(waveFreqs),length(win));
    for bc = 1:length(waveFreqs)
        curData = waveDecomp(bc,win);
        pow(bc,:) = real(curData).^2 + imag(curData).^2;
        mpow(bc) = 10*log10(mean(pow(bc,:)));
        phaseList(bc,:) = atan(imag(curData)./real(curData));
    end
    powersHigh{ck} = pow;
    mPowerHigh{ck} = mpow;
    phasesHigh{ck} = phaseList;
end

%Now low powers 
for ck = 1:length(lowIDXList)
    curLoc = lowIDXList(ck);
    if curLoc > round(fs*PowerCalculationWindow/2)
        if curLoc<length(asymmetry)-round(fs*PowerCalculationWindow/2)
            win = curLoc-round(fs*PowerCalculationWindow/2):curLoc+round(fs*PowerCalculationWindow/2)-1;
        else
            win = curLoc-round(fs*PowerCalculationWindow/2)-(length(asymmetry)-curLoc-1):length(asymmetry);
        end
    else
        win = 1:curLoc+(round(fs*PowerCalculationWindow/2))+(round(fs*PowerCalculationWindow/2)-curLoc);
    end
%     if curLoc > round(fs*PowerCalculationWindow/2)
%         if curLoc<length(asymmetry)-round(fs*PowerCalculationWindow/2)
%             win = curLoc-round(fs*PowerCalculationWindow/2):curLoc+round(fs*PowerCalculationWindow/2)-1;
%         else
%             win = curLoc-round(fs*PowerCalculationWindow/2)-(length(asymmetry)-curLoc-1):length(asymmetry);
%         end
%     else
%         win = 1:curLoc+(round(fs*PowerCalculationWindow/2)-1)+(curLoc-round(fs*PowerCalculationWindow/2)-1);
%     end
    pow = zeros(length(waveFreqs),length(win));
    mpow = zeros(1,length(waveFreqs));
    phaseList = zeros(length(waveFreqs),length(win));
    for bc = 1:length(waveFreqs)
        curData = waveDecomp(bc,win);
        pow(bc,:) = real(curData).^2 + imag(curData).^2;
        mpow(bc) = 10*log10(mean(pow(bc,:)));
        phaseList(bc,:) = atan(imag(curData)./real(curData));
    end
    powersLow{ck} = pow;
    mPowerLow{ck} = mpow;
    phasesLow{ck} = phaseList;
end
