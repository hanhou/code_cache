%% Extract raw spikes from imec .ap.bin file guided by kilosort times
% (for Karel 20210823) 
% Using https://djoshea.github.io/neuropixel-utils/analysis/

%% Set paths
rawFile = 'H:\HH08\HH08_S02_20210813_g0\HH08_S02_20210813_g0_imec0\HH08_S02_20210813_g0_t0.imec0.ap.bin';
ks2Folder = 'I:\catGT\HH08\catgt_HH08_S02_20210813_g0\HH08_S02_20210813_g0_imec0\imec0_ks2';

chanMap = 'D:\Han_Sync\Svoboda\Scripts\Ephys\neuropixel-utils\map_files\neuropixPhase3B1_kilosortChanMap.mat';

%% Load data
addpath('D:\Han_Sync\Svoboda\Scripts\Ephys\neuropixel-utils');

imec = Neuropixel.ImecDataset(rawFile, 'channelMap', chanMap);
ks = Neuropixel.KilosortDataset(ks2Folder, 'imecDataset', imec, 'channelMap', chanMap);
ks.load();

%% Drift map
ks.printBasicStats();
stats = ks.computeBasicStats()
metrics = ks.computeMetrics();
metrics.plotDriftmap('driftThreshold', 4, 'spikeAmpQuantile', 0.8);

%% Cluster drift map
cluster_ids = metrics.cluster_ids(1:5:metrics.nClusters);
metrics.plotClusterDriftmap('cluster_ids', cluster_ids, 'colorByAmp', true, ...
                            'showSmooth', true, 'smoothWidthSeconds', 500, ...
                            'showIndividual', true, 'onlyGaps', true, ...
                            'zShuffleClusters', true, ...
                            'alpha', 0.9, 'scaleWithAmp', false);

%% Plot cluster centers of mass
metrics.plotClusterWaveformAtCentroid()

%% Raw waveforms
% https://djoshea.github.io/neuropixel-utils/waveforms/
cluster_id = 100;

snippetSet = ks.getWaveformsFromRawData('cluster_ids', cluster_id, 'window', [-30, 60], ...
    'num_waveforms', 100, 'best_n_channels', 20, 'car', true, ... 
    'subtractOtherClusters', true, 'centerUsingFirstSamples', 30);
snippetSet.plotAtProbeLocations('alpha', 0.5);

%% Export some good spike snippsets for Karel
nCell = 20;
nWaveform = 100;
window = [-20, 40];

% Get best nCell cells
[~, bestInds] = sort(metrics.cluster_amplitude, 'descend');
bestInds = intersect(bestInds, ks.clusters_good, 'stable');  % Only keep good KS units
figure();

allSpikesuV = nan(range(window) + 1, nWaveform, nCell);  % [nTime, nLength, nCell]
for i = 1:nCell
    ax = subplot(4,5,i);
    disp(i);
    cluster_id = metrics.cluster_ids(bestInds(i));
    snippetSet = ks.getWaveformsFromRawData('cluster_ids', cluster_id, 'window', window, ...
                    'num_waveforms', nWaveform, 'best_n_channels', 1, 'car', true, ... 
                    'subtractOtherClusters', true, 'centerUsingFirstSamples', 20); 
    allSpikesuV(:, :, i) = double(squeeze(snippetSet.data)) * snippetSet.scaleToUv(1);
    plot(ax, snippetSet.time_ms', allSpikesuV(:,:,i), 'k');
    drawnow; 
end

time_ms = snippetSet.time_ms;

%% Save data
toSave = 1:20 ; % [1, 2, 6, 7, 9, 10, 12, 13, 15, 19, 20];
spikes = allSpikesuV(:,:,toSave);
save('spike_snippets.mat', 'spikes', 'time_ms');








