load('psthToDave.mat');
figure(); 

% Artifact removal
artifactEdges = [0, 0.002];  % 0 and 2 ms
artifactSpan = [-0.125, 0.3] * 1e-3; % Relative to the eaah edge
artifacts = false(size(allSpikeTimeAligned));
for i = 1:length(artifactSpan)
    thisArtifacts = artifactEdges(i)+artifactSpan(1) <= allSpikeTimeAligned & allSpikeTimeAligned <= artifactEdges(i)+artifactSpan(2);
    artifacts = artifacts | thisArtifacts;
end
allSpikeTimeAligned(artifacts) = [];
allSpikeDepth(artifacts) = [];

% Fig 8D
figDDepthROI = [1.5 2.5];   % Only CC in figure D
figDTimeBinSize = 0.2 * 1e-3; % 0.25 ms
figDTimeBins = -0.01:figDTimeBinSize:0.02;   % 0.25 ms win
figDSpikeTimes = allSpikeTimeAligned(figDDepthROI(1) <= allSpikeDepth & allSpikeDepth <= figDDepthROI(2));
[n, ~] = hist(figDSpikeTimes, figDTimeBins);

psthEvoked = n/nReps/figDTimeBinSize;
psthBaseline = mean(psthEvoked(figDTimeBins < 0));

subplot(1,2,1); 
plot(figDTimeBins*1000, psthEvoked - psthBaseline, 'k-');

hold on; plot([0 0], ylim, 'k--'); plot([2 2], ylim, 'k--'); plot(xlim, [0 0], 'k--')
xlabel('Time from stimulus onset (ms)');
ylabel('delta spikes/s');
title(sprintf('time bin = %g ms', figDTimeBinSize*1000));

% Fig 8E
figEDepthBins = 0:0.02:4;  % 20 um win
figETimeROI = [0 0.02];   % Evoked spikes
figEBaselineROI = [-0.01 0];   % Baseline

figEDepths = allSpikeDepth(figETimeROI(1) <= allSpikeTimeAligned & allSpikeTimeAligned <= figETimeROI(2));
[n, ~] = hist(figEDepths, figEDepthBins); 
psthEvoked = n / nReps / range(figETimeROI);

figEDepthsBaseline = allSpikeDepth(figEBaselineROI(1) <= allSpikeTimeAligned & allSpikeTimeAligned <= figEBaselineROI(2));
[n, ~] = hist(figEDepthsBaseline, figEDepthBins); 
psthBaseline = n / nReps / range(figEBaselineROI);

subplot(1,2,2); 
plot(psthEvoked - psthBaseline, figEDepthBins, 'k-');
title(sprintf('depth bin = %g um', (figEDepthBins(2) - figEDepthBins(1))*1000));
hold on; plot([0 0], ylim, 'k--')
set(gca,'YDir','reverse')
xlabel('delta spikes/s');
ylabel('Depth in the brain (mm)');

