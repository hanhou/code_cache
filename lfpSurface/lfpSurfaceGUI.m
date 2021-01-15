function [lfpCorrCoef, lfpSurfaceChAuto] = lfpSurfaceGUI(lfpFilename, lfpPara)
%% Estimate the brian surface from LFP channels, using 10 segments of 1 sec each samples. 
% Part of code is based on https://github.com/cortex-lab/spikes
% Han Hou 2021

nClips = 10;
clipDurMax = 1; % seconds
startTime = 3; % skip first seconds

lfpFs = lfpPara.Fs;  % sampling rate of lfp
nChansInFile = lfpPara.nChansInFile; 
antiStaggering = lfpPara.antiStaggering;  % Average LFP from the sites of the same row on the probe
freqBandForSurface = lfpPara.freqBandForSurface;  % Freq band for lfp power
outChannels = lfpPara.outChannels; % Channels that are outside the brain for sure

% Locate nClips samples
d = dir(lfpFilename); 
nSamps = d.bytes/2/nChansInFile;
nClipSamps = round(lfpFs*clipDurMax);
sampStarts = round(linspace(lfpFs*startTime, nSamps-nClipSamps, nClips+1)); 
nClipSamps = min(round(lfpFs*clipDurMax), diff(sampStarts(1:2)));  % Ensure the clips are not over-lapped (for short file)

mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});

chanIdx = lfpPara.chanIdx;  % Meaningful channels on one shank

if antiStaggering
    nDepth = length(chanIdx)/2; % Anti-stagerring in LFP due to the layout of NPX probes
else
    nDepth = length(chanIdx);
end

rawDat = zeros(nDepth, nClipSamps, nClips);
lfpCorrCoef = zeros(nDepth, nDepth, nClips);
lfpCov = zeros(nDepth, nDepth, nClips);
powerForSurface = zeros(nDepth, nClips);

% For each clips
for n = 1:nClips
%     disp(n)
    
    % pull out the data
    thisDat = double(mmf.Data.x(chanIdx, (1:nClipSamps)+sampStarts(n)));  % Get data and remove the sync channel
    
    if ~isempty(lfpPara.refChan)
        thisDat(lfpPara.refChan,:) = thisDat(lfpPara.refChan-1, :); % Deal with reference channel
    end
    
    if lfpPara.lowPassFilter
        parfor ch = 1:size(thisDat,1)
            thisDat(ch,:) = lowpass(thisDat(ch,:), 300, lfpFs, 'ImpulseResponse','iir','Steepness',0.8);
        end
    end

    % median subtract 
    thisDat = bsxfun(@minus, thisDat, median(thisDat,2));  % Self-median
    if ~isempty(outChannels)
        thisDat = thisDat - mean(thisDat(outChannels, :),1);   % Median of channels that are outside the brain for sure (good for removing 60 Hz)
    end
    
    % Spatial averaging across the same depths (two sites in a row)
    if antiStaggering
        thisDat = squeeze(mean(reshape(thisDat',size(thisDat,2),2,[]),2))';  
    end
        
    % Power spectrum
    [Pxx, F] = myTimePowerSpectrumMat(thisDat', lfpFs);
    if n==1
        allPowerEst = zeros(size(Pxx,1), size(Pxx,2), nClips);
    end
    allPowerEst(:,:,n) = Pxx;

    % Low band power for LFP surface
    inclF = F>freqBandForSurface(1) & F<=freqBandForSurface(2);
    powerForSurface(:, n) = 10*log10(mean(Pxx(inclF,:)));
    
    % Correlation
    lfpCorrCoef(:, :, n) = corrcoef(thisDat');
    lfpCov(:, :, n) = cov(thisDat');
    
    % Cache raw data
    rawDat(:,:,n) = thisDat;
end

% -- Average over clips --
lfpCorrAver = median(lfpCorrCoef, 3);
lfpCovAver = median(lfpCov, 3);
powerForSurfaceAver = median(powerForSurface, 2);
allPowerEst = mean(allPowerEst, 3);
[lfpSurfaceChAuto, ~] = findSurface(powerForSurfaceAver, lfpCorrAver);
if antiStaggering
    lfpSurfaceChAuto = lfpSurfaceChAuto * 2;
end

% -- Generate plots --
fig = figure(); hold on;
set(gcf,'defaultAxesFontSize',15);
set(gcf,'uni','norm','pos',[0.019       0.054       0.957       0.867]);
ax(1) = subplot(2,3,[1 4]);
ax(2) = subplot(2,3,2);
ax(3) = subplot(2,3,3);
ax(4) = subplot(2,3,5);
ax(5) = subplot(2,3,6);

% -- Attach data --
myData.lfpFilename = lfpFilename;
myData.lfpPara = lfpPara;
myData.lfpSurfaceCh = lfpSurfaceChAuto;
myData.F = F;
myData.allPowerEst = allPowerEst;
myData.powerForSurface = powerForSurfaceAver;
myData.lfpCorr = lfpCorrAver;
myData.lfpCov = lfpCovAver;
myData.rawDat = rawDat;
myData.ax = ax;
myData.ylimUser = [];

% -- Save button --
myData.saveBtn = uicontrol('Style', 'pushbutton', 'String', 'Save user-defined surface channel', ...
    'Callback', @saveBtn , ...
    'Fontsize',15, 'uni', 'norm', 'pos', [0.004 0.942 0.15 0.047]);
set(myData.saveBtn, 'Enable', 'off');

% -- Load saved user surface, if any --
[path, ~, ~] = fileparts(lfpFilename);
if exist(fullfile(path, 'lfpSurfaceUser.txt'))
    myData.lfpSurfaceChUser = dlmread(fullfile(path, 'lfpSurfaceUser.txt'));
else
    myData.lfpSurfaceChUser = nan; 
end

set(fig, 'UserData', myData);
updateFig(fig);


function [lfpSurfaceChAuto, lfpCorrForSurface] = findSurface(lfpPower, lfpCorrCoef)
% Use power & correlation to automatically (but roughly) find the surface channel
corrAverRange = 30;  % Range to average the corr matrix, from the initial guess
initialGuessByPower = find(lfpPower > median(lfpPower), 1, 'last'); % Last channel with power > median (in the brain for sure)
corrToAver = initialGuessByPower-corrAverRange: initialGuessByPower;
lfpCorrForSurface = mean(lfpCorrCoef(corrToAver, :));  % Average the correlation coeff
[~, lfpSurfaceChAuto] = min(diff(smooth(lfpCorrForSurface, 10)));  % Fastest decay of the corr coeff

function [Pxx, F] = myTimePowerSpectrumMat(x, Fs)
L = size(x,1);
NFFT = 2^nextpow2(L);
[Pxx,F] = pwelch(x,[],[],NFFT,Fs);

function updateFig(fig)
myData = get(fig, 'UserData');
lfpCorr = myData.lfpCorr;
lfpCov = myData.lfpCov;
depthOnProbe = (0 : size(lfpCorr)-1)' * (1 + myData.lfpPara.antiStaggering) * myData.lfpPara.siteDistance/2; % in um

% -- LFP raw around the surface --
drawRaw(depthOnProbe, myData.lfpPara.Fs, myData.rawDat, myData.lfpSurfaceCh, myData.lfpSurfaceChUser, myData.ax(1));

% -- Power --
drawPower(depthOnProbe, myData.F, myData.allPowerEst, myData.lfpSurfaceCh, myData.lfpSurfaceChUser, myData.ax([2, 3]));

% -- Corr coef and Covariance --
drawCorrCov(depthOnProbe, lfpCorr, myData.lfpSurfaceCh, myData.lfpSurfaceChUser, myData.ax(4))
title(myData.ax(4), 'LFP correlation');
drawCorrCov(depthOnProbe, lfpCov, myData.lfpSurfaceCh, myData.lfpSurfaceChUser, myData.ax(5))
title(myData.ax(5), 'LFP covariance');

for a = 2:5
    set(zoom(myData.ax(a)), 'ActionPostCallback', @updateZoom);
end

drawnow;
linkaxes(myData.ax(2:end), 'y')
linkaxes([myData.ax(4) myData.ax(5)], 'x')

function drawRaw(depthOnProbe, samplingRate, rawDat, lfpSurfaceCh, lfpSurfaceChUser, ax, gain)
myData = get(gcf, 'UserData');
nChBelow = 20;
nChAbove = 10;
if ~isnan(lfpSurfaceChUser)
    [~, surfaceIdx] = min(abs(depthOnProbe - lfpSurfaceChUser * myData.lfpPara.siteDistance/2));
else
    [~, surfaceIdx] = min(abs(depthOnProbe - lfpSurfaceCh * myData.lfpPara.siteDistance/2));
end
drawRange = max(1, surfaceIdx-nChBelow): min(surfaceIdx+nChAbove, size(rawDat,1));

% -- Reorganize data --
rawToDraw = rawDat(drawRange, :, :);
rawToDraw = reshape(rawToDraw, size(rawToDraw,1), [], 1);

axes(ax); cla; hold on;
if nargin < 7
    gain = 1/range(rawToDraw(:)) * mean(diff(depthOnProbe)) * 3;
end

for dd = 1:length(drawRange)
    offset = depthOnProbe(drawRange(dd));
    thisLFP = rawToDraw(dd,:); 
    thisLFP = thisLFP * gain;
    plot((1:length(thisLFP))*1000/samplingRate, -thisLFP + offset, 'b');
end
xlabel('Time (ms)')
axis tight; ylimTmp = ylim();
drawSurfaceLines(lfpSurfaceCh, lfpSurfaceChUser);
xlim([0 2000]);
ylim(ylimTmp);
ylabel('Distance to probe tip (um)');

% -- Annotation --
title(sprintf('LFP surface channel:\nauto = %g, \\color{red}user-defined = %g\n', ...
    myData.lfpSurfaceCh, lfpSurfaceChUser));

% Add callback
set(ax, 'ButtonDownFcn', {@newUserSurface});

function drawCorrCov(depthOnProbe, data, lfpSurfaceCh, lfpSurfaceChUser, ax)
myData = get(gcf, 'UserData');
axes(ax); cla;
h_corr = imagesc(depthOnProbe, depthOnProbe, data); hold on;
colorbar; colormap jet;
set(gca,'YDir','normal') 
title('LFP covariance')
xlabel('Distance to probe tip (um)');
ylabel('Distance to probe tip (um)');

% Surface lines
h_l = drawSurfaceLines(lfpSurfaceCh, lfpSurfaceChUser, 2);

if ~isempty(myData.ylimUser); ylim(myData.ylimUser);  end

% Add callback
set([h_corr; h_l(:)], 'ButtonDownFcn', {@newUserSurface});

function drawPower(depthOnProbe, F, allPowerEst, lfpSurfaceCh, lfpSurfaceChUser, axs)
myData = get(gcf, 'UserData');

dispRange = [0 100]; % Hz
freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};
dispF = F>dispRange(1) & F<=dispRange(2);

% -- Power spectrum --
axes(axs(1)); cla
h_power = imagesc(F(dispF), depthOnProbe, 10*log10(allPowerEst(dispF,:))');
xlim(dispRange);
xlabel('Frequency (Hz)');
ylabel('Distance to probe tip (um)');
set(gca, 'YDir', 'normal'); hold on;
title('Power specdtrum')
h_l1 = drawSurfaceLines(lfpSurfaceCh, lfpSurfaceChUser);
if ~isempty(myData.ylimUser); ylim(myData.ylimUser);  end

% -- 1-d Power --
axes(axs(2)); cla; hold on;
c = copper(length(freqBands));
c = c(:, [3 2 1]);
set(axs(2), 'ColorOrder', c);
for q = 1:length(freqBands)
    inclF = F>freqBands{q}(1) & F<=freqBands{q}(2);
    thisPow = mean(10*log10(allPowerEst(inclF,:)),1);
    plot(thisPow, depthOnProbe, 'linewidth', 2);
end
h_l2 = drawSurfaceLines(lfpSurfaceCh, lfpSurfaceChUser);
legend([cellfun(@(x)sprintf('%.1f - %.1f Hz', x(1), x(2)), freqBands, 'uni', false), 'auto', 'user-defined'],...
    'Location', 'southwest', 'Fontsize', 15);
xlabel('Power (dB)');
ylabel('Distance to probe tip (um)');
ylim([min(depthOnProbe) max(depthOnProbe)])
title('Power per freq band')
if ~isempty(myData.ylimUser); ylim(myData.ylimUser);  end

% Add callback
set([h_power; h_l1(:); h_l2(:); axs(2)], 'ButtonDownFcn', {@newUserSurface});

function h = drawSurfaceLines(lfpSurfaceCh, lfpSurfaceChUser, ~)
myData = get(gcf, 'UserData');
if nargin>2; drawX = true; else; drawX = false; end
% Surface lines
posTmp = [lfpSurfaceCh lfpSurfaceCh] * myData.lfpPara.siteDistance/2;
h(1) = plot(xlim(), posTmp, 'k--', 'linew', 2);
if drawX; h(2) = plot(posTmp, ylim(), 'k--', 'linew', 2); end
if ~isnan(lfpSurfaceChUser)
    posTmp = [lfpSurfaceChUser lfpSurfaceChUser]*myData.lfpPara.siteDistance/2;
    plot(xlim(), posTmp, 'r', 'linew', 2);
    if drawX; plot(posTmp, ylim(), 'r', 'linew', 2); end
end

function newUserSurface(~, ~)
% Update user-defined surface channel
pos = get(gca,'currentPoint');
y = pos(1,2);
myData = get(gcf, 'UserData');
myData.lfpSurfaceChUser = round(y/(myData.lfpPara.siteDistance/2));
set(gcf, 'UserData', myData);
updateFig(gcf);
set(myData.saveBtn, 'Enable', 'on');

function updateZoom(~, ~)
myData = get(gcf, 'UserData');
myData.ylimUser = ylim();

% Sync all axis limits
for a=2:5
    set(myData.ax(a), 'ylim', myData.ylimUser);
    if a>=4
        set(myData.ax(a), 'xlim', myData.ylimUser);
    end
end
set(gcf, 'UserData', myData);

function saveBtn(~, ~)
myData = get(gcf, 'UserData');
[path, ~, ~] = fileparts(myData.lfpFilename);
if ~isnan(myData.lfpSurfaceChUser)
    dlmwrite(fullfile(path, 'lfpSurfaceUser.txt'), myData.lfpSurfaceChUser);
end
set(myData.saveBtn, 'Enable', 'off');

