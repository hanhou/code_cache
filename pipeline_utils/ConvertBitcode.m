% Based on Dave's code
function ConvertBitcode(sessionDir)
%% Settings
if nargin == 0
    sessionDir=uigetdir('F:\catGT\');
end

% In line with my pybpod code
bitCodeDigits = 20;

chan.protocol = 2;  % Protocol dig marker channel
eventMarkerDur.bitcodeEachbit = 1;
eventMarkerDur.goCue = 10;
eventMarkerDur.reward = 20;
eventMarkerDur.ITI = 30;

chan.behavior = 2;  % Behavior dig marker channel (I moved all behavior channel to bitcode channel in Apr 21 version)
eventMarkerDur.choiceL = 2;   % I made a mistake in Apr version that choiceL and the start of bitcode are the SAME....
eventMarkerDur.choiceR = 3;

%% get Ephys Bitcode
% Start of a trial (onset of my bitcode is indicated by twice of the bitcode width)
sTrig = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.bitcodeEachbit * 2)); 

% Solve the bitcode-start and choiceL overlapping bug...
ITI = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.ITI));  
choiceL = [];
for i=1:length(ITI)
    this_choice_idx = find(sTrig(i) < sTrig & sTrig < ITI(i));
    choiceL = [choiceL; sTrig(this_choice_idx)];
    sTrig(this_choice_idx) = [];
end

eTrig=[sTrig(2:end); sTrig(end)+20];
bitsAll = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.bitcodeEachbit)); 
goCue = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.goCue)); 
reward = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.reward)); 
% choiceL = dlmread(getTxtFileName(sessionDir, chan.behavior, eventMarkerDur.choiceL)); 
choiceR = dlmread(getTxtFileName(sessionDir, chan.behavior, eventMarkerDur.choiceR)); 

bitcode = zeros(length(sTrig), bitCodeDigits);

for i = 1:length(sTrig)  % For each trial
    bitHighThis = bitsAll(bitsAll>sTrig(i) & bitsAll<eTrig(i));
    bitHighPositionThis = round((bitHighThis - sTrig(i) - eventMarkerDur.bitcodeEachbit * 1e-3) / (2 * eventMarkerDur.bitcodeEachbit * 1e-3));
    bitcode(i, bitHighPositionThis) = 1;
end

bitCodeS = num2str(bitcode, '%d');

%% Save .mat files for ingestion
imecFolders = dir(fullfile(sessionDir, '*imec*'));
for f = 1:length(imecFolders)  % Save the same bitcode.mat to each imec folder 
    toSave = [imecFolders(f).name(1 : strfind(imecFolders(f).name,'imec')-1) 'bitcode.mat'];
    save(fullfile(imecFolders(f).folder, imecFolders(f).name, toSave), ...
        'bitcode', 'bitCodeS', 'goCue', 'sTrig', 'reward', 'choiceL', 'choiceR', 'ITI');
end

function txtFile = getTxtFileName(sessionDir, chan, duration)
txtFile = dir(fullfile(sessionDir, sprintf('*%g_%g.adj.txt', chan, duration)));  % Try TPrime adjusted first
if isempty(txtFile)  % If no adj.txt, try raw txt
    txtFile = dir(fullfile(sessionDir, sprintf('*%g_%g.txt', chan, duration)));
end    
txtFile = fullfile(txtFile.folder, txtFile.name);
    
