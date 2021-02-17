

function [event_trig_lfp, event_trig_lfp_aver, event_trig_CSD, ts] = pstLFPByDepth(lfpFilename, lfpFs, nChansInFile, eventTimes, lfpSurfaceCh, antiStaggering)
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
    thisLFP = double(mmf.Data.x(1 : end-1, startIdx : min(end, startIdx+length(ts)-1)));  % Remove sync channel
    
    thisLFP(192,:) = mean(thisLFP(191, :),1); % Duplicate reference channel
    thisLFP = bsxfun(@minus, thisLFP, median(thisLFP,2));  % Remove self-median
    thisLFP = thisLFP - mean(thisLFP(lfpSurfaceCh+3:end, :),1);   % Remove median of air channels (remove 60 Hz)

    if antiStaggering
        thisLFP = squeeze(mean(reshape(thisLFP',size(thisLFP,2),2,[]),2))';  % Spatial averaging across the same depths (two sites in a row)
%         thisLFP = thisLFP(1:2:end, :, :);  % Only use half of sites
    end
    
    event_trig_lfp(:,1:size(thisLFP,2),ee) = thisLFP;
end

% Postprocessing
gainFactor = 2.3438; 
event_trig_lfp = event_trig_lfp * gainFactor;  % In uV

% -- Compute CSD --
event_trig_lfp_aver = nanmedian(event_trig_lfp, 3);
event_trig_lfp_aver = smoothdata(event_trig_lfp_aver,1, 'movmedian', 5);

event_trig_CSD = CSD(event_trig_lfp_aver' * 1e-6, lfpFs, (1+antiStaggering)*10e-6, 'ifplot', 0, 'inverse', 2*(1+antiStaggering)*10e-6)';
event_trig_CSD = smoothdata(event_trig_CSD, 1, 'gaussian', 5);
% event_trig_CSD = CSD(event_trig_lfp_aver' * 1e-6, lfpFs, (1+antiStaggering)*10e-6, 'ifplot', 1)';

