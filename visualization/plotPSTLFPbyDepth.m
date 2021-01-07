
function plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType, ax, lfpSurfaceCh, chNumToShow)
% function plotPSTHbyDepth(timeBins, depthBins, allP, eventName, ax)
%
% see matching function "psthByDepth" that generates input for this
%
% psthType is 'norm' if you have normalized units, otherwise 

nD = length(depthBins)-1;

if nargin<6 || isempty(ax)    
    ax = gca;
end

% imagesc(timeBins, depthBins(1:end-1), allP);
actualDepth = lfpSurfaceCh * 10 - depthBins(1:end-1); 
imagesc(timeBins, actualDepth, allP);
set(gca, 'YDir', 'normal');
hold on;

if strcmp(psthType, 'norm')
%     plot(zeros(1, nD), depthBins(1:end-1), 'k--', 'LineWidth', 2.0)
    colormap(ax, colormap_BlueWhiteRed)
    caxis([-10 10]);
else
%     plot(zeros(1, nD), depthBins(1:end-1), 'w--', 'LineWidth', 2.0)
end
    
xlim([-500 1400]);
ylim([-500 3900]);

xlabel(['time from ' eventName ' (ms)']);
% xlim([min(timeBins) max(timeBins)]);
% ylabel('depth on electrode array (µm)')
ylabel('Distance from pia (LFP surface, µm)')
box off
h = colorbar;
if strcmp(psthType, 'norm')
    h.Label.String = 'Firing rate z-score';
else
    h.Label.String = 'Firing rate (spikes/sec)';
end

% Show channel number as texts
for i = 1:length(chNumToShow)
    text(min(xlim), (lfpSurfaceCh - chNumToShow(i)) * 10, sprintf('ch %g', chNumToShow(i)), 'color', 'k')
end

% makepretty;
set(gca, 'YDir','reverse')
