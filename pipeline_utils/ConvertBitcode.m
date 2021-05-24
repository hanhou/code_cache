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

bitcode = zeros(length(iti), bitCodeDigits);   % Use length iti to make sure (the last) trial has ended.

[STRIG_, GOCUE_, CHOICEL_, CHOICER_, REWARD_, ITI_] = deal(1,2,3,4,5,6);

% Use iti as trial marker to exclude truncated trials
digMarkerPerTrial = nan(length(iti), 6);  % [STrig, goCue, choiceL, choiceR, reward, ITI]

% digMarkerPerTrial(:,[STRIG_, GOCUE_, ITI_]) = [sTrig, goCue, iti];  % Must exists 
digMarkerPerTrial(:,ITI_) = iti;  % Must exists 

for i = 1:length(iti)  % For each trial
    thisStart = sTrig(find(sTrig < iti(i), 1, 'last'));  % Deal with truncated trials (only search backward from each iti)
    digMarkerPerTrial(i, STRIG_) = thisStart;
    digMarkerPerTrial(i, GOCUE_) = goCue(thisStart < goCue & goCue < iti(i));
    
    % Parse bitcode
    bitHighThis = bitsAll(thisStart < bitsAll & bitsAll < iti(i));
    bitHighPositionThis = round((bitHighThis - thisStart - eventMarkerDur.bitcodeEachbit * 1e-3) / (2 * eventMarkerDur.bitcodeEachbit * 1e-3));
    
    % Error-proof
    assert(max(bitHighPositionThis) <= bitCodeDigits, 'bitCodeDecodingError');
    
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

% Reassign goCue and STrig to make sure they're all complete trials (i.e., paired with an iti)
sTrig = digMarkerPerTrial(:,STRIG_);
goCue = digMarkerPerTrial(:,GOCUE_);

% Debug
% figure(); plot(goCue,1,'g>'); hold on; plot(iti,2,'k*'); plot(sTrig,0,'bo'); ylim([-2 5]);
% for i = 1:length(iti)    ;text(iti(i),2,num2str(i));  end
% for i = 1:length(goCue)   ;text(goCue(i),1,num2str(i)); end
% for i = 1:length(sTrig)   ;text(sTrig(i),0,num2str(i)); end


%% Override trial num fix (in the case where len(ephys) > len(behav))
trialNumFixName = {'HH09_20210406', 'HH09_20210418'};
trialNumFix = {1:164, 1:513};

%% Save .mat files for ingestion
imecFolders = dir(fullfile(sessionDir, '*imec*'));
for f = 1:length(imecFolders)  % Save the same bitcode.mat to each imec folder 
    % For DJ pipeline
    matName = [imecFolders(f).name(1 : strfind(imecFolders(f).name,'imec')-1) 'bitcode.mat'];
    fullFileNameDJ = fullfile(imecFolders(f).folder, imecFolders(f).name, matName);
    
    % Trial num fix?
    fixIndex = find(~cellfun(@isempty, regexp(sessionDir, trialNumFixName)));
    if ~isempty(fixIndex)
        trialNum = trialNumFix{fixIndex};
        bitcode = bitcode(trialNum);
        bitCodeS = bitCodeS(trialNum,:);
        goCue = goCue(trialNum);
        sTrig = sTrig(trialNum);
        digMarkerPerTrial = digMarkerPerTrial(trialNum,:);
        save(fullFileNameDJ, ...
            'bitcode', 'bitCodeS', 'goCue', 'sTrig', 'reward', 'choiceL', 'choiceR', 'iti', 'digMarkerPerTrial', 'trialNum');
        fprintf('Trial Number fixed!!\n')
    else
        save(fullFileNameDJ, ...
            'bitcode', 'bitCodeS', 'goCue', 'sTrig', 'reward', 'choiceL', 'choiceR', 'iti', 'digMarkerPerTrial');
    end
    
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
    
