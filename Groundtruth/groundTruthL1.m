%% Ground truth experiment analysis
clear
% -- Settings --
% MU thresholding
MUThr = 0;
gainFactor = 0.6/512/500*1e6;
antiLFPStaggering = true;

% Depth PSTH
depthBinSize = 20; % in units of the channel coordinates, in this case µm
timeBinSize = 0.001; % seconds
bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization

ksDir = uigetdir('e:\catGT\');

% cat_folder = 'E:\catGT\HH100\catgt_GroundTruth04_g0\';
% run_name = 'GroundTruth04';
% probe = 'imec0';
% cat_folder = 'E:\catGT\HH101\catgt_GroundTruth06_dual_bilateral_g0\';
% run_name = 'GroundTruth06_dual_bilateral';
% probe = 'imec0';

%%
% -- Fetch digital signals --
% trialStart = dlmread([cat_folder run_name '_g0_tcat.nidq.XD_2_0_0.adj.txt']); % Left CC
% goCue = dlmread([cat_folder run_name '_g0_tcat.nidq.XD_2_1_1200.adj.txt']); % Left CC

% -- Fetch Kilosort results--
% ksDir = [cat_folder run_name '_g0_' probe '\' probe '_ks2\'];

%{
amp=readNPY([ks_folder 'amplitudes.npy']);
% spktime=double(readNPY([ks_folder 'spike_times_sec.npy']));
spktime=double(readNPY([ks_folder 'spike_times.npy']))/30000;
channel_positions = readNPY([ks_folder 'channel_positions.npy']);

% spike_clusters = readNPY([ks_folder 'spike_clusters.npy']);
% temps = readNPY([ks_folder 'templates.npy']);
% winv = readNPY([ks_folder 'whitening_mat_inv.npy']);  % Should be winv!!!
% ycoords = channel_positions(:,2);
%}

trialStartFile = dir(fullfile(ksDir, '..\..\', '*XD_2_0_0.adj.txt'));
goCueFile = dir(fullfile(ksDir, '..\..\', '*XD_2_1_1200.adj.txt'));

trialStart = dlmread(fullfile(trialStartFile.folder, trialStartFile.name));
goCue = dlmread(fullfile(goCueFile.folder, goCueFile.name));
trialN = length(goCue); % trial number of the spike

% Compensate drift of goCue VS photostim (due to catGT error) !!!
if exist(fullfile(ksDir, 'catGTLastError.txt'))
    catGTLastError = dlmread(fullfile(ksDir, 'catGTLastError.txt'));
    fprintf('------------- catGT ERROR corrected: %g sec/%g trials ---------\n', catGTLastError, trialN)
else
    catGTLastError = 0; 
end
goCue = goCue + (1:trialN)' * catGTLastError/trialN;
photoStimTime = goCue + 1.2;

% Use loadKSdir
sp = loadKSdir(ksDir);

% Get spike depth and amplitude in uV
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
spikeAmpsuV = spikeAmps * gainFactor;

%{
cluN=readNPY('H:\groundtruthL1\HH101\3\spike_sites.npy');
thr=40;
cluN=cluN(amp>thr);
bigU=unique(cluN);
spktime=spktime(amp>thr);
trialN=ones(size(spktime))*length(goCue); % trial number of the spike
spktime2=spktime;
%% get the spiketimes (new JRClust) for photostim latency
cluN=double(centerSites(:,2)); % cluster number is the site number
bigU=unique(cluN);
spktime=double(spikeTimes)/30000;
trialN=ones(size(spktime))*length(goCue); % trial number of the spike
spktime2=spktime;
%}

% Thresholding
bigUnits = spikeAmpsuV >= MUThr;
sp.st = sp.st(bigUnits);
sp.clu = sp.clu(bigUnits);
spikeDepths = spikeDepths(bigUnits);

spktime = sp.st;
uniqueMUPosition = unique(spikeDepths);
spktime2 = spktime;

lfpSurfaceCh = 384;  % Dummy

%% --- Draw LFP, get surface channel from LFP ---
lfpFile = dir(fullfile(ksDir, '..\', '*.lf.bin'));
lfpFileFullname = fullfile(lfpFile.folder, lfpFile.name);

lfpFs = 2500;  % neuropixels phase3a
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

% -- Compute LFP power and find surface channel --
[lfpByChannel, allPowerEst, F, allPowerVar, lfpCorr, lfpSurfaceCh] = ...
    lfpBandPower(lfpFileFullname, lfpFs, nChansInFile, [], [0, 20], antiLFPStaggering);

% chanMap = readNPY(fullfile(ksDir, 'channel_map.npy'));
% nC = length(chanMap);
nC = nChansInFile - 1;
depthOnProbe = (0:nC-1)*10;

if antiLFPStaggering
    nC = nC/2;
end
allPowerEst = allPowerEst';

% plot LFP power and LFP surface
dispRange = [0 100]; % Hz
marginalChans = round(linspace(10, nC, 6));
freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};

plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands, lfpSurfaceCh, antiLFPStaggering);
set(gcf, 'name', ksDir)
set(gcf,'uni','norm','pos',[0.338       0.145       0.443       0.693]);
SetFigure(20);

% plot LFP correlation
figure('Name', 'LFP correlation'); imagesc(depthOnProbe, depthOnProbe, lfpCorr); hold on;
plot(xlim(), [lfpSurfaceCh lfpSurfaceCh]*10, 'r--', 'linew', 2);
plot([lfpSurfaceCh lfpSurfaceCh]*10, ylim(), 'r--', 'linew', 2);
colorbar; colormap gray;
set(gca,'YDir','normal') 
SetFigure(20);
title(sprintf('LFP surface = %g um', lfpSurfaceCh*10));

%% --- stim-triggered LFP ---
photoStimTime = photoStimTime + 0.0014; % LFP offset

fprintf('Stim-triggered LFP.....');
[event_trig_lfp_all, event_trig_lfp_aver, event_trig_CSD, lfp_ts] = pstLFPByDepth(lfpFileFullname, lfpFs, nChansInFile, photoStimTime, lfpSurfaceCh, antiLFPStaggering);
fprintf('Done!\n');

%% -- Plotting
figure('name', ksDir)
set(gcf,'uni','norm','pos',[0.009        0.08        0.64       0.826]);

ax1_lfp = subplot(1,2,1); 
plotLFPbyDepth(lfp_ts * 1000, event_trig_lfp_aver, ax1_lfp, lfpSurfaceCh, [1, 384], antiLFPStaggering)
hold on; 
plot(xlim(), [0 0], 'b--', 'linew', 2);
text(min(xlim())+100, -100, sprintf('LFP surface: ch #%g', lfpSurfaceCh), 'color', 'b');
plot(xlim(), lfpSurfaceCh*10 - [3840 3840], 'k-');
plot(xlim(), [lfpSurfaceCh*10 lfpSurfaceCh*10], 'k-');

otherTimeMarkers = [0:200:1000 (0:200:1000)+2];
for i = 1:length(otherTimeMarkers)
    x = otherTimeMarkers(i);
    plot(x*[1 1], ylim(), 'k--');
end
title(sprintf('Evoked LFP, %g trials', trialN));

ax2_lfp = subplot(1,2,2); 
plotLFPbyDepth(lfp_ts * 1000, event_trig_CSD, ax1_lfp, lfpSurfaceCh, [1, 384], antiLFPStaggering)
hold on; 
%colormap(ax2_lfp, jet)
linkaxes([ax1_lfp, ax2_lfp], 'xy');
plot(xlim(), [0 0], 'b--', 'linew', 2);
text(min(xlim())+100, -100, sprintf('LFP surface: ch #%g', lfpSurfaceCh), 'color', 'b');
plot(xlim(), lfpSurfaceCh*10 - [3840 3840], 'k-');
plot(xlim(), [lfpSurfaceCh*10 lfpSurfaceCh*10], 'k-');

for i = 1:length(otherTimeMarkers)
    x = otherTimeMarkers(i);
    plot(x*[1 1], ylim(), 'k--');
end
title(sprintf('Evoked CSD, %g trials', trialN));
% colormap(gca, flipud(jet)); % blue = source; red = sink
title('CSD, \muA/mm^3 (\color{blue}sink, \color{red}source\color{black})');

SetFigure(20);
drawnow;

% -- Show some averaged LFP
figure('name', ['Stim-triggered LFP', ksDir]);
set(gcf,'uni','norm','pos',[0.651       0.081       0.339       0.825]);
hold on;
lfpActualDepthToShow = -300:300:3500;
lfpActualDepth = (lfpSurfaceCh - (0:nC-1)*(1+antiLFPStaggering))*10; 
gain = 5;

for dd = 1:length(lfpActualDepthToShow)
    % subplot(length(lfpActualDepthToShow), 1, dd);
    [~, thisIdx] = min(abs(lfpActualDepthToShow(dd) - lfpActualDepth));
    offset = lfpActualDepth(thisIdx);
    thisLFP = squeeze(event_trig_lfp_all(thisIdx,:,:))';
    thisLFP = thisLFP / range(thisLFP(:)) * mean(diff(lfpActualDepthToShow) * gain);
    shadedErrorBar(lfp_ts * 1000, -thisLFP + offset, {@median, @(x) 2*std(x,[],1)/sqrt(size(x,1))});
end

for i = 1:length(otherTimeMarkers)
    x = otherTimeMarkers(i);
    plot(x*[1 1], ylim(), 'k--');
end

% linkaxes(findall(gcf,'type','axes'),'y');
linkaxes([gca ax1_lfp],'xy');
set(gca, 'YDir','reverse')
title(sprintf('%g trials', trialN));
xlim([-200 1400])
SetFigure(20);

%% --- psthViewer --
% psthViewer(spikeTimes, clu, eventTimes, window, trGroups)
% if your events come in different types, like different orientations of a
% visual stimulus, then you can provide those values as "trial groups",
% which will be used to construct a tuning curve. Here we just give a
% vector of all ones. 
trialGroups = ones(size(sp.events)); 

window = [0, 3];
otherTimeMarkers = [1.2:0.2:2.2 (1.2:0.2:2.2)+0.002];
psthViewer(sp.st, sp.clu, goCue, window, trialGroups, templateDepths, waveforms, otherTimeMarkers, lfpSurfaceCh)
set(gcf, 'name', ksDir)


%% --- PSTHbyDepth ---
psthType = ''; % show the normalized version
eventName = 'stimulus onset'; % for figure labeling

[timeBins, depthBins, allP, normVals] = psthByDepth(sp.st, spikeDepths, ...
    depthBinSize, timeBinSize, goCue, window, bslWin);
timeBins = (timeBins - 1.2) * 1000;

figure('name', ksDir)
set(gcf,'uni','norm','pos',[0.142       0.073        0.78       0.825]);
ax1_spk = subplot(1,2,1); 
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, 'norm', [], lfpSurfaceCh, [1, 384]);  % Normalized

hold on; 
plot(xlim(), [0 0], 'b--', 'linew', 2);
text(min(xlim())+100, -100, sprintf('LFP surface: ch #%g', lfpSurfaceCh), 'color', 'b');
plot(xlim(), lfpSurfaceCh*10 - [3840 3840], 'k-');
plot(xlim(), [lfpSurfaceCh*10 lfpSurfaceCh*10], 'k-');
title(sprintf('Evoked spikes, %g trials', trialN));

otherTimeMarkers = [0:200:1000 (0:200:1000)+2];
for i = 1:length(otherTimeMarkers)
    x = otherTimeMarkers(i);
    plot(x*[1 1], ylim(), 'k--');
end

ax2_spk = subplot(1,2,2); 
plotPSTHbyDepth(timeBins, depthBins, allP .* normVals(:,2), eventName, '', [], lfpSurfaceCh, [1, 384]);   % Only mean is removed
hold on; 
colormap(ax2_spk, hot)
linkaxes([ax1_lfp, ax2_lfp, ax1_spk, ax2_spk], 'xy');
plot(xlim(), [0 0], 'b--', 'linew', 2);
text(min(xlim())+100, -100, sprintf('LFP surface: ch #%g', lfpSurfaceCh), 'color', 'b');
plot(xlim(), lfpSurfaceCh*10 - [3840 3840], 'k-');
plot(xlim(), [lfpSurfaceCh*10 lfpSurfaceCh*10], 'k-');

for i = 1:length(otherTimeMarkers)
    x = otherTimeMarkers(i);
    plot(x*[1 1], ylim(), 'k--');
end
SetFigure(20);
drawnow;

%% --- Driftmap ---
[spikeTimes, spikeAmps, spikeDepths_driftmap, spikeSites] = ksDriftmap(ksDir);
figure('name', ksDir)
plotDriftmap(spikeTimes, spikeAmps, spikeDepths_driftmap);
plot(xlim(), [lfpSurfaceCh lfpSurfaceCh]*10, 'b--', 'linew', 2);
text(min(xlim())+100, lfpSurfaceCh*10+100, sprintf('LFP surface: ch #%g', lfpSurfaceCh), 'color', 'b');

% makepretty;

%% --- Spike Amps ---
depthBins = 0:40:3840;
ampBins = 0:30:min(max(spikeAmps),800);
recordingDur = sp.st(end);

[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths_driftmap, ampBins, depthBins, recordingDur);
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins, lfpSurfaceCh);

%% --- Retrieve any unit from the plot ---
yPos = 1725;
ind = abs(uniqueMUPosition - yPos) < depthBinSize;
closest_pos = uniqueMUPosition(ind);
for i=1:length(closest_pos)
    fprintf('Closest #%g: cid = %g, depth = %g\n', i, mode(sp.clu(spikeDepths == closest_pos(i))), closest_pos(i));
end
fprintf('\n')


