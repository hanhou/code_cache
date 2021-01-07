

function [event_trig_lfp, ts] = pstLFPByDepth(lfpFilename, lfpFs, nChansInFile, eventTimes, lfpSurfaceCh, antiStaggering)
% function [timeBins, depthBins, allP] = psthByDepth(spikeTimes, ...
%   spikeDepths, depthBinSize, timeBinSize, eventTimes, win, bslWin[, bslEvents])
%
% Computes LFP triggered by events and current source density
 
%% load LFP raw data
d = dir(lfpFilename); 
nSamps = d.bytes/2/nChansInFile;
mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});

nChan= nChansInFile - 1; 
if antiStaggering
    nDepth = nChan/2;
else
    nDepth = nChan;
end

nTrial = length(eventTimes);
win = [-0.5, 2];  % Window in sec
ts = win(1): 1/lfpFs: win(2);
event_trig_lfp = nan(nDepth, length(ts), nTrial);

for ee = 1 : nTrial
    startIdx = round((eventTimes(ee)+win(1))*lfpFs);
    thisLFP = double(mmf.Data.x(1 : end-1, startIdx : startIdx+length(ts)-1));  % Remove sync channel
    
    thisLFP(192,:) = mean(thisLFP(191, :),1); % Duplicate reference channel
    thisLFP = bsxfun(@minus, thisLFP, median(thisLFP,2));  % Remove self-median
    thisLFP = thisLFP - mean(thisLFP(lfpSurfaceCh+3:end, :),1);   % Remove median of air channels (remove 60 Hz)

    if antiStaggering
        thisLFP = squeeze(mean(reshape(thisLFP',size(thisLFP,2),2,[]),2))';  % Spatial averaging across the same depths (two sites in a row)
    end
    
    event_trig_lfp(:,:,ee) = thisLFP;
end
