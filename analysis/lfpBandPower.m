

function [lfpByChannel, allPowerEst, F, allPowerVar, lfpCorr, lfpSurfaceCh] = lfpBandPower(lfpFilename, lfpFs, nChansInFile, freqBands, freqBandForSurface)
% function [lfpByChannel, allPowerEst, F] = lfpBandPower(lfpFilename, lfpFs, nChansInFile, freqBand)
% Computes the power in particular bands, and across all frequencies, across the recording
% samples 10 segments of 1 sec each to compute these things. 

if ~isempty(freqBands) && ~iscell(freqBands)
    freqBands = {freqBands};
end
nF = length(freqBands);

nClips = 10;
clipDur = 1; % seconds
startTime = 3; % skip first seconds
% surfacePowerCutoff = 3; % Surface channel = Drop -3 dB power from median (50% median power)
corrAverRange = 10;

% load nClips one-sec samples
d = dir(lfpFilename); 
nSamps = d.bytes/2/nChansInFile;
sampStarts = round(linspace(lfpFs*startTime, nSamps, nClips+1)); % skip first 10 secs
nClipSamps = round(lfpFs*clipDur);

mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});

lfpRaw = mmf.Data.x(1:end-1, :);  % Remove sync channel
nChansInFile = nChansInFile - 1;
lfpRaw(192,:) = lfpRaw(191, :); % Duplicate reference channel

allPowerEstByBand = zeros(nClips, nChansInFile, nF);
lfpCovs = zeros(nChansInFile, nChansInFile, nClips);
powerLowBand = zeros(nChansInFile, nClips);

for n = 1:nClips
%     fprintf(1, 'clip%d: ', n);
    % pull out the data
    thisDat = double(lfpRaw(:, (1:nClipSamps)+sampStarts(n)));
    
    % median subtract? 
    thisDat = bsxfun(@minus, thisDat, median(thisDat,2));  % Self-median
%     thisDat = bsxfun(@minus, thisDat, mean(thisDat,2));
    thisDat = thisDat - mean(thisDat(380:384, :),1);   % Median of air channels (remove 60 Hz)
        
    [Pxx, F] = myTimePowerSpectrumMat(thisDat', lfpFs);
    
    if n==1
        allPowerEst = zeros(nClips, size(Pxx,1), size(Pxx,2));
    end
    allPowerEst(n,:,:) = Pxx;
        
    for f = 1:nF
        inclF = F>freqBands{f}(1) & F<=freqBands{f}(2);
        allPowerEstByBand(n,:, f) = mean(Pxx(inclF,:));
    end
    
    % Cache power ([0,20]) across channels
    inclF = F>freqBandForSurface(1) & F<=freqBandForSurface(2);
    powerLowBand(:, n) = 10*log10(mean(Pxx(inclF,:)));
%     brainPowerLevel = median(powerForSurface); % Median of all channels (assuming >50% channels are in the brain)
%     surfaceChannels(n) = find(powerForSurface < brainPowerLevel - surfacePowerThreshold, 1);
%     fprintf('surface ch = %g\n', surfaceChannels(n));
%     plot(powerForSurface, 'k');
%     plot([surfaceChannels(n) surfaceChannels(n)], ylim(), 'r--')
    
    % Cache LFP correlation
    lfpCovs(:, :, n) = corrcoef(thisDat');
    
end

if nF>0
    lfpByChannel = squeeze(mean(allPowerEstByBand, 1)); % mean across clips
else
    lfpByChannel = [];
end
allPowerVar = squeeze(var(allPowerEst,1));
allPowerEst = squeeze(mean(allPowerEst, 1));

lfpCorr = mean(lfpCovs, 3);

% -- Use power & correlation to find the surface channel --
powerLowBand = median(powerLowBand, 2);
surfaceGuessByPower = find(powerLowBand > median(powerLowBand), 1, 'last'); % Last channel with power > median (in the brain for sure)
corrToAver = surfaceGuessByPower-corrAverRange: surfaceGuessByPower;
lfpCorrForSurface = mean(lfpCorr(corrToAver, :));  % Average the correlation coeff
[~, lfpSurfaceCh] = min(diff(smooth(lfpCorrForSurface)));  % Fastest decay of the corr coeff



function [Pxx, F] = myTimePowerSpectrumMat(x, Fs)
L = size(x,1);
NFFT = 2^nextpow2(L);
[Pxx,F] = pwelch(x,[],[],NFFT,Fs);