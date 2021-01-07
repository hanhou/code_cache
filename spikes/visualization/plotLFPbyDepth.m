
function plotLFPbyDepth(timeBins, lfp, ax, lfpSurfaceCh, chNumToShow, antiStaggering)
% function plotPSTLFPbyDepth(timeBins, depthBins, allP, eventName, ax)
%
% see matching function "pstLFPByDepth" that generates input for this
%

if nargin<6 || isempty(ax)    
    ax = gca;
end

% imagesc(timeBins, depthBins(1:end-1), allP);

depthOnProbe = (0:size(lfp,1)-1)*10;
if antiStaggering
    depthOnProbe = depthOnProbe * 2;
end

actualDepth = lfpSurfaceCh * 10 - depthOnProbe; 

imagesc(timeBins, actualDepth, lfp);
set(gca, 'YDir','reverse')

xmax = max(lfp(:)); xmin = min(lfp(:));
colormap_BlueWhiteRed(300, 1.5, -(xmax + xmin)/(xmax-xmin));
hold on;
colorbar;
    
xlim([-500 1400]);
ylim([-500 3900]);

xlabel(['time (ms)']);
ylabel('Distance from pia (LFP surface, µm)')
box off

% Show channel number as texts
for i = 1:length(chNumToShow)
    text(min(xlim), (lfpSurfaceCh - chNumToShow(i)) * 10, sprintf('ch %g', chNumToShow(i)), 'color', 'k')
end

%makepretty;

