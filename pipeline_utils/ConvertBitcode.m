% Based on Dave's code
function ConvertBitcode(sessionDir)
%% Settings
if nargin == 0
    sessionDir=uigetdir('F:\catGT\');
end

% In line with my pybpod code
bitCodeDigits = 20;

chan.protocol = 1;  % Protocol dig marker channel
eventMarkerDur.bitcodeEachbit = 1;
eventMarkerDur.goCue = 10;
eventMarkerDur.reward = 20;
eventMarkerDur.ITI = 0;

chan.behavior = 0;  % Behavior dig marker channel (I moved all behavior channel to bitcode channel in Apr 21 version)
eventMarkerDur.choiceL = 1;   % I made a mistake in Apr version that choiceL and the start of bitcode are the SAME....
eventMarkerDur.choiceR = 2;

%% get Ephys Bitcode
% Start of a trial (onset of my bitcode is indicated by twice of the bitcode width)
sTrig = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.bitcodeEachbit * 2)); 

% Solve the bitcode-start and choiceL overlapping bug...
iti_all = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.ITI));  
% choiceL = [];
% for i=1:length(iti)
%     this_choice_idx = find(sTrig(i) < sTrig & sTrig < iti(i));
%     choiceL = [choiceL; sTrig(this_choice_idx)];
%     sTrig(this_choice_idx) = [];
% end

bitsAll = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.bitcodeEachbit)); 
goCue = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.goCue)); 
reward = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.reward)); 
choiceL = dlmread(getTxtFileName(sessionDir, chan.behavior, eventMarkerDur.choiceL)); 
choiceR = dlmread(getTxtFileName(sessionDir, chan.behavior, eventMarkerDur.choiceR)); 

iti = setdiff(iti_all, [bitsAll; goCue; reward; sTrig]);

bitcode = zeros(length(sTrig), bitCodeDigits);

[STRIG_, GOCUE_, CHOICEL_, CHOICER_, REWARD_, ITI_] = deal(1,2,3,4,5,6);

% Use iti as trial marker to exclude truncated trials
digMarkerPerTrial = nan(length(iti), 6);  % [STrig, goCue, choiceL, choiceR, reward, ITI]

% digMarkerPerTrial(:,[STRIG_, GOCUE_, ITI_]) = [sTrig, goCue, iti];  % Must exists 
digMarkerPerTrial(:,ITI_) = iti;  % Must exists 

for i = 1:length(iti)  % For each trial
    thisStart = sTrig(find(sTrig < iti(i), 1, 'last'));  % Deal with truncated trials
    digMarkerPerTrial(i, STRIG_) = thisStart;
    digMarkerPerTrial(i, GOCUE_) = goCue(thisStart < goCue & goCue < iti(i));
    
    % Parse bitcode
    bitHighThis = bitsAll(thisStart < bitsAll & bitsAll < iti(i));
    bitHighPositionThis = round((bitHighThis - sTrig(i) - eventMarkerDur.bitcodeEachbit * 1e-3) / (2 * eventMarkerDur.bitcodeEachbit * 1e-3));
    bitcode(i, bitHighPositionThis) = 1;
    
    % Fill in the digMarkerPerTrial matrix
    thisL = choiceL(thisStart < choiceL & choiceL < iti(i));
    if ~isempty(thisL)
        digMarkerPerTrial(i, CHOICEL_) = thisL;
    end
    thisR = choiceR(thisStart < choiceR & choiceR < iti(i));
    if ~isempty(thisR)
        digMarkerPerTrial(i, CHOICER_) = thisR;
    end        
    thisRew = reward(thisStart < reward & reward < iti(i));
    if ~isempty(thisRew)
        digMarkerPerTrial(i, REWARD_) = thisRew;
    end             
end

bitCodeS = num2str(bitcode, '%d');

%% Save .mat files for ingestion
imecFolders = dir(fullfile(sessionDir, '*imec*'));
for f = 1:length(imecFolders)  % Save the same bitcode.mat to each imec folder 
    % For DJ pipeline
    matName = [imecFolders(f).name(1 : strfind(imecFolders(f).name,'imec')-1) 'bitcode.mat'];
    fullFileNameDJ = fullfile(imecFolders(f).folder, imecFolders(f).name, matName);
    save(fullFileNameDJ, ...
        'bitcode', 'bitCodeS', 'goCue', 'sTrig', 'reward', 'choiceL', 'choiceR', 'iti', 'digMarkerPerTrial');
    
    % For phy event_plugin_hh
    ks2Folder = dir(fullfile(fullfile(imecFolders(f).folder, imecFolders(f).name), 'imec*_ks2'));
    fullFileNamePhy = fullfile(ks2Folder.folder, ks2Folder.name, matName);
    copyfile(fullFileNameDJ, fullFileNamePhy);
    
    fprintf('%s saved!\n', fullFileNameDJ);
    fprintf('%s saved!\n', fullFileNamePhy);
    
end

function txtFile = getTxtFileName(sessionDir, chan, duration)
txtFile = dir(fullfile(sessionDir, sprintf('*%g_%g.adj.txt', chan, duration)));  % Try TPrime adjusted first
if isempty(txtFile)  % If no adj.txt, try raw txt
    txtFile = dir(fullfile(sessionDir, sprintf('*%g_%g.txt', chan, duration)));
    fprintf('No _adj.txt found, using non-adjusted version');
end    
txtFile = fullfile(txtFile.folder, txtFile.name);
    
