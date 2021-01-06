%% Ground truth experiment analysis
clear
% -- Settings --
% MU thresholding
MUThr = 0;
gainFactor = 0.6/512/500*1e6;

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
trialN = ones(size(spktime))*length(goCue); % trial number of the spike
spktime2 = spktime;

%% realign the spiketimes for each trial
viT_offset_file2=[trialStart' trialStart(end)+20]-0.5;
for i = length(viT_offset_file2)-1:-1:1
    trialN(spktime>=viT_offset_file2(i) & spktime<viT_offset_file2(i+1))=i;
    spktime2(spktime>=viT_offset_file2(i) & spktime<viT_offset_file2(i+1))=(spktime(spktime>=viT_offset_file2(i) & spktime<viT_offset_file2(i+1))-goCue(i));
end

% re-zero
spktime2=spktime2-1.2; % for the pulse stimulation (only sample cue)

% get psth to photostim, Plot latency to photo-stim
% close all
trialMax=max(trialN);

binwidth=0.002; % s
bins=-3:binwidth:3; % bins

ignoreFR=zeros(length(bins)-1,length(uniqueMUPosition)); % for CD
stimFR=zeros(length(bins)-1,length(uniqueMUPosition)); % for CD
timePeriod=[-1.2 -1.19]; % Before trial start
timePeriodS=[0 .1]; % get a function of time period and depth
diffFR=zeros(1,length(uniqueMUPosition)); % for CD

for i=1:length(uniqueMUPosition)
    nSpk= histcounts(spktime2(spikeDepths==uniqueMUPosition(i)),bins)'/length(trialStart)/binwidth;    
    nSpkS=histcounts(spktime2(spikeDepths==uniqueMUPosition(i)),bins)'/length(trialStart)/binwidth;
    
    ignoreFR(:,i)=nSpk; % for CD
    stimFR(:,i)=nSpkS; % for CD
    
    sSpk=histcounts(spktime2(spikeDepths==uniqueMUPosition(i)),timePeriod)'/length(trialStart)/binwidth;    
    sSpkS=histcounts(spktime2(spikeDepths==uniqueMUPosition(i)),timePeriodS)'/length(trialStart)/binwidth;
    
    diffFR(i)=sSpkS-sSpk;
end

%% plot
%
% close all
xmax=100;
xmin=-30;
% % plot(diffFR,(length(diffFR):-1:1)*3.84/length(diffFR))
% figure; 
% plot(diffFR,double(uniqueMUPosition))
% box off
% title('0-20ms')
% ylabel('Probe length (mm)')
% xlabel('Difference in MUA (spike/s)')
% set(gca,'TickDir','out')
% ax=gca;
% % ax.YDir = 'reverse';
% %ylim([0 3.84]);

%%
figure('name', ksDir)
% timePeriod=[-1.23 -1.21];
timePeriod=[-0.04 -0.02]; % pulses
startI=find(round(bins,2)==round(timePeriod(1),2));
endI=find(round(bins,2)==round(timePeriod(2),2));
% diffFR2=mean(stimFR(startI:endI,:),1)-mean(ignoreFR(startI:endI,:),1);
diffFR2=mean(stimFR(startI:endI,:),1)-mean(stimFR(1:500,:),1); % pulses
subplot(131)
plot(diffFR2,double(uniqueMUPosition))
box off
title('-20-0 ms')
% ylabel('Depth in the brain (mm)')
ylabel('Length on the probe (mm)')
ax=gca;
% ax.YDir = 'reverse';
%ylim([0 3.84])
xlim([xmin xmax])

% timePeriod=[-1.21 -1.19];
stims=zeros(6,length(uniqueMUPosition));
for i = 1:6
    timePeriod=[0.006 0.026]+(i-1)*0.2; % pulses
    startI=find(round(bins,2)==round(timePeriod(1),2));
    endI=find(round(bins,2)==round(timePeriod(2),2));
    stims(i,:)=mean(stimFR(startI:endI,:),1)-mean(stimFR(1:500,:),1);
end

% diffFR2=mean(stimFR(startI:endI,:),1)-mean(ignoreFR(startI:endI,:),1);
diffFR20=mean(stims);
subplot(132)
xE=double(uniqueMUPosition);
plot(diffFR20,xE)
box off
title('0-20 ms')
set(gca,'TickDir','out')
ax=gca;
% ax.YDir = 'reverse';
%ylim([0 3.84])
xlim([xmin xmax])
xlabel('Difference in MUA (spike/s)')

% timePeriod=[-1.19 -1.17];
timePeriod=[.02 .04]; % pulses
startI=find(round(bins,2)==round(timePeriod(1),2));
endI=find(round(bins,2)==round(timePeriod(2),2));
% diffFR2=mean(stimFR(startI:endI,:),1)-mean(ignoreFR(startI:endI,:),1);
diffFR2=mean(stimFR(startI:endI,:),1)-mean(stimFR(1:500,:),1);
subplot(133)
plot(diffFR2,double(uniqueMUPosition))
box off
title('20-40 ms')
set(gca,'TickDir','out')
ax=gca;
% ax.YDir = 'reverse';
xlim([xmin xmax])
%ylim([0 3.84])

%% --- psthViewer --
% psthViewer(spikeTimes, clu, eventTimes, window, trGroups)
% if your events come in different types, like different orientations of a
% visual stimulus, then you can provide those values as "trial groups",
% which will be used to construct a tuning curve. Here we just give a
% vector of all ones. 
trialGroups = ones(size(sp.events)); 

window = [0, 3];
otherTimeMarkers = [1.2:0.2:2.2 (1.2:0.2:2.2)+0.002];
psthViewer(sp.st, sp.clu, sp.events, window, trialGroups, otherTimeMarkers)

%% --- PSTHbyDepth ---
psthType = ''; % show the normalized version
eventName = 'stimulus onset'; % for figure labeling

[timeBins, depthBins, allP, normVals] = psthByDepth(sp.st, spikeDepths, ...
    depthBinSize, timeBinSize, sp.events, window, bslWin);
timeBins = (timeBins - 1.2) * 1000;

figure('name', ksDir)
set(gcf,'uni','norm','pos',[0.142       0.073        0.78       0.825]);
ax1 = subplot(1,2,1); 
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, 'norm');   % Normalized
ylim([0 3500])
hold on; 
otherTimeMarkers = [0:200:1000 (0:200:1000)+2];
for i = 1:length(otherTimeMarkers)
    x = otherTimeMarkers(i);
    plot(x*[1 1], ylim(), 'k--');
end

ax2 = subplot(1,2,2); 
plotPSTHbyDepth(timeBins, depthBins, allP .* normVals(:,2), eventName, '');   % Only mean is removed
ylim([0 3500])
hold on; 
for i = 1:length(otherTimeMarkers)
    x = otherTimeMarkers(i);
    plot(x*[1 1], ylim(), 'w--');
end
colormap(ax2, jet)
linkaxes([ax1, ax2], 'xy');


%% --- Driftmap ---
[spikeTimes, spikeAmps, spikeDepths_driftmap, spikeSites] = ksDriftmap(ksDir);
figure('name', ksDir)
plotDriftmap(spikeTimes, spikeAmps, spikeDepths_driftmap);

%% --- Spike Amps ---
depthBins = 0:40:3840;
ampBins = 0:30:min(max(spikeAmps),800);
recordingDur = sp.st(end);

[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths_driftmap, ampBins, depthBins, recordingDur);
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);

%% --- LFP ---
%{
lfDir = [cat_folder run_name '_g0_' probe '\'];
lfpD = dir(fullfile(lfDir, '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(lfDir, lfpD(1).name);

lfpFs = 2500;  % neuropixels phase3a
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

[lfpByChannel, allPowerEst, F, allPowerVar] = ...
    lfpBandPower(lfpFilename, lfpFs, nChansInFile, []);

chanMap = readNPY(fullfile(ksDir, 'channel_map.npy'));
nC = length(chanMap);

allPowerEst = allPowerEst(:,chanMap+1)'; % now nChans x nFreq

% plot LFP power
dispRange = [0 100]; % Hz
marginalChans = [10:50:nC];
freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};

plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands);

%}
%% --- Retrieve any unit from the plot ---
yPos = 2580;
[~,id] = min(abs(uniqueMUPosition - yPos));
cid = mode(sp.spikeTemplates(spikeDepths == uniqueMUPosition(id)))


