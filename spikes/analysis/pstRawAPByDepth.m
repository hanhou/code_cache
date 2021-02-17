

function [event_trig_raw, event_trig_raw_aver, ts] = pstRawAPByDepth(rawFilename, rawFs, nChansInFile, eventTimes)
% function [timeBins, depthBins, allP] = psthByDepth(spikeTimes, ...
%   spikeDepths, depthBinSize, timeBinSize, eventTimes, win, bslWin[, bslEvents])
%
% Computes LFP triggered by events and current source density
 
%% load LFP raw data
d = dir(rawFilename); 
nSamps = d.bytes/2/nChansInFile;
mmf = memmapfile(rawFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});

nChan= nChansInFile - 1; 
nTrial = length(eventTimes);
win = [-0.01, 0.05];  % Window in sec
ts = win(1): 1/rawFs: win(2);
event_trig_raw = nan(nChan, length(ts), nTrial);

for ee = 1 : nTrial
    startIdx = round((eventTimes(ee)+win(1))*rawFs);
    thisRaw = double(mmf.Data.x(1 : end-1, startIdx : startIdx+length(ts)-1));  % Remove sync channel
    thisRaw = bsxfun(@minus, thisRaw, mean(thisRaw,2));  % -<T>

    event_trig_raw(:,:,ee) = thisRaw;
end

% Postprocessing
gainFactor = 2.3438; 
event_trig_raw = event_trig_raw * gainFactor;  % In uV
event_trig_raw_aver = median(event_trig_raw, 3);

