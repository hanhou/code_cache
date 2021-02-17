% Estimate the brian surface from LFP channels, using 10 segments of 1 sec each samples. 
% Part of code is based on https://github.com/cortex-lab/spikes
% Han Hou 2021

clear 
% close all

%% == Parameters ==
probe = '1.0';
% probe = '2.0';

if strcmp(probe, '1.0')
    lfpPara.Fs = 2500;  % sampling rate of lfp
    lfpPara.lowPassFilter = false;  % Additional low pass filter
    lfpPara.outChannels = 380:384;  % Channels that are outside the brain for sure
    lfpPara.nChansInFile = 385;   % Including the sync channel, if any  
    lfpPara.chanIdx = 1:384;   % All channels on one shank 
    lfpPara.refChan = 192;
    lfpPara.siteDistance = 20; % vertical distance of the recording sites (um)
    lfpPara.antiStaggering = true;  % Average LFP from the sites of the same row on the probe
    lfpPara.freqBandForSurface = [0.5 20];  % Freq band for lfp power
else  % NPX2.0 data (not fully tested)
    lfpPara.Fs = 30000;  % sampling rate of lfp
    lfpPara.lowPassFilter = true;  % Additional low pass filter
    lfpPara.nChansInFile = 385;   % Including the sync channel, if any  
    lfpPara.siteDistance = 15; % vertical distance of the recording sites (um)
    lfpPara.chanIdx = [1:48, 97:144];   % Shank 02 on the four-shank probe
    lfpPara.refChan = [];
    lfpPara.antiStaggering = false;  % Average LFP from the sites of the same row on the probe
    lfpPara.freqBandForSurface = [0.5 20];  % Freq band for lfp power
    lfpPara.outChannels = [];  % Channels that are outside the brain for sure
end

%% -- Get file --
[lfpFile, lfpPath] = uigetfile('*.bin');
lfpFileFullname = fullfile(lfpPath, lfpFile);

%% -- Plot lfp surface --
[lfpCorr, lfpSurfaceCh_auto] = lfpSurfaceGUI(lfpFileFullname, lfpPara);