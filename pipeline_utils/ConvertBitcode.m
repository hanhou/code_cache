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
eventMarkerDur.bitcodeFirstMultiplier = 1.5;  % The multiplier of the first bit code width
eventMarkerDur.goCue = 10;
eventMarkerDur.reward = 20;
eventMarkerDur.ITI = 30;

chan.behavior = 2;  % Behavior dig marker channel (I moved all behavior channel to bitcode channel in Apr 21 version)
eventMarkerDur.choiceL = 2;  
eventMarkerDur.choiceR = 2.5;

chan.bpodstart = 1;
chan.zaber = 4;
chan.cameras = [5, 6, 7];

%% get Ephys Bitcode
% Start of a trial (onset of my bitcode is indicated by twice of the bitcode width)
sTrig = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.bitcodeEachbit * eventMarkerDur.bitcodeFirstMultiplier)); 
iti = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.ITI));

bitsAll = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.bitcodeEachbit));
goCue = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.goCue));
reward = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.reward));
choiceL = dlmread(getTxtFileName(sessionDir, chan.behavior, eventMarkerDur.choiceL));
choiceR = dlmread(getTxtFileName(sessionDir, chan.behavior, eventMarkerDur.choiceR));
choiceR = setdiff(choiceR, choiceL);   % 2.5 * (1-20%) = 2.0!!! choice R also includes all choice L...

bpodTrialStart = dlmread(getTxtFileName(sessionDir, chan.bpodstart, 0));
zaberStepsUp = dlmread(getTxtFileName(sessionDir, chan.zaber, 0));
zaberStepsDwn = dlmread(getTxtFileName(sessionDir, chan.zaber, 0, '.iXD'));
zaberStepsAll = sort([zaberStepsUp; zaberStepsDwn]);

for i = 1:length(chan.cameras)
    cameras{i} = dlmread(getTxtFileName(sessionDir, chan.cameras(i), 0));
end

bitcode = zeros(length(iti), bitCodeDigits);   % Use length iti to make sure (the last) trial has ended.

[STRIG_, GOCUE_, CHOICEL_, CHOICER_, REWARD_, ITI_, BPOD_START_, ZABER_IN_POS_] = deal(1,2,3,4,5,6,7,8);

% Use iti as trial marker to exclude truncated trials
digMarkerPerTrial = nan(length(iti), 8);  % [STrig, goCue, choiceL, choiceR, reward, ITI]
zaberPerTrial = {};
cameraPerTrial = {};

% digMarkerPerTrial(:,[STRIG_, GOCUE_, ITI_]) = [sTrig, goCue, iti];  % Must exists 
digMarkerPerTrial(:,ITI_) = iti;  % Must exists 

for i = 1:length(iti)  % Fill in digMarkerPerTrial
    thisBitCodeStart = sTrig(find(sTrig < iti(i), 1, 'last'));  % Deal with truncated trials (only search backward from each iti)
    digMarkerPerTrial(i, BPOD_START_) = bpodTrialStart(find(bpodTrialStart < iti(i), 1, 'last'));   % The closest bpodTrialStart before this ITI
    digMarkerPerTrial(i, STRIG_) = thisBitCodeStart;
    digMarkerPerTrial(i, GOCUE_) = goCue(thisBitCodeStart < goCue & goCue < iti(i));
        
    % Parse bitcode
    bitHighThis = bitsAll(thisBitCodeStart < bitsAll & bitsAll < iti(i));
    bitHighPositionThis = round((bitHighThis - thisBitCodeStart - (eventMarkerDur.bitcodeFirstMultiplier - 1) * eventMarkerDur.bitcodeEachbit * 1e-3) ...
                                 / (2 * eventMarkerDur.bitcodeEachbit * 1e-3));
    
    % Bitcode bugfix
    if contains(sessionDir, 'HH09_S06_20210609') && i == 97
        bitHighPositionThis = [1     4     5     7    10    13    15    17    20];
    elseif contains(sessionDir, 'HH09_S08_20210612') && i == 250
        bitHighPositionThis = [ 1     3     4     6     7    11    13    14    15    16    18    19    20];
    elseif contains(sessionDir, 'HH09_S08_20210612') && i == 329
        bitHighPositionThis = [1     3     7     8     9    10    11    12    14    18];
    elseif contains(sessionDir, 'HH09_S09_20210613') && i == 495
        bitHighPositionThis = [1     2     4     6     8     9    10    11    13    15    17    19];
    else
        assert(max(bitHighPositionThis) <= bitCodeDigits, 'bitCodeDecodingError');
    end

   
    % Set bitcode
    bitcode(i, bitHighPositionThis) = 1;
    
    % Fill in the digMarkerPerTrial matrix
    thisL = choiceL(thisBitCodeStart < choiceL & choiceL < iti(i));
    if ~isempty(thisL)
        digMarkerPerTrial(i, CHOICEL_) = thisL;
    end
    thisR = choiceR(thisBitCodeStart < choiceR & choiceR < iti(i));
    if ~isempty(thisR)
        digMarkerPerTrial(i, CHOICER_) = thisR;
    end        
    thisRew = reward(thisBitCodeStart < reward & reward < iti(i));
    if ~isempty(thisRew)
        digMarkerPerTrial(i, REWARD_) = thisRew;
    end
    
    % zaber steps (only forward protraction pulses)
    zaberPerTrial{i} = zaberStepsAll(thisBitCodeStart < zaberStepsAll & zaberStepsAll < iti(i));
    digMarkerPerTrial(i, ZABER_IN_POS_) = max(zaberPerTrial{i});   % The last zaber pulse of this trial
     
end

for i = 1:length(iti)  % Fill in camera pulses
    thisStart = digMarkerPerTrial(i, BPOD_START_);  % Use BPOD_START instead of STrig to capture camera pulses during ITI
    if i < length(iti)
        thisEnd = digMarkerPerTrial(i + 1, BPOD_START_);
    else
        thisEnd = Inf;
    end
       
    for cc = 1:length(chan.cameras)
        cameraPerTrial{cc}{i} = cameras{cc}(thisStart < cameras{cc} & cameras{cc} < thisEnd);
    end
end

bitCodeS = num2str(bitcode, '%d');

% Reassign goCue and STrig to make sure they're all complete trials (i.e., paired with an iti)
sTrig = digMarkerPerTrial(:,STRIG_);
goCue = digMarkerPerTrial(:,GOCUE_);
bpodTrialStart = digMarkerPerTrial(:,BPOD_START_); 

choiceAll = digMarkerPerTrial(:, [CHOICEL_, CHOICER_]);
choiceAll = choiceAll(:);
choiceAll = sort(choiceAll(~isnan(choiceAll)));

% % Debug
% figure(); plot(goCue,1,'g>'); hold on; plot(iti,2,'k*'); plot(sTrig,0,'bo'); ylim([-2 5]);
% plot(zaberPerTrial{1},4,'k+'); plot(cameraPerTrial{1},3,'k.');
% for i = 1:length(iti)    ;text(iti(i),2,num2str(i));  end
% for i = 1:length(goCue)   ;text(goCue(i),1,num2str(i)); end
% for i = 1:length(sTrig)   ;text(sTrig(i),0,num2str(i)); end


%% Override trial num fix (in the case where len(ephys) > len(behav))
trialNumFixName = {'HH09_20210406', 'HH09_20210418', 'HH09_20210429'};
trialNumFix = {1:164, 1:513, 1:393};

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
            'bitcode', 'bitCodeS', 'goCue', 'sTrig', 'reward', 'choiceL', 'choiceR', 'choiceAll', 'iti', 'digMarkerPerTrial', ...
            'bpodTrialStart', 'zaberStepsUp', 'zaberStepsDwn', 'cameras', ...
            'chan', 'eventMarkerDur', ...
            'zaberPerTrial', 'cameraPerTrial', ...
            'trialNum');
        fprintf('Trial Number fixed!!\n')
    else
        save(fullFileNameDJ, ...
            'bitcode', 'bitCodeS', 'goCue', 'sTrig', 'reward', 'choiceL', 'choiceR', 'choiceAll', 'iti', 'digMarkerPerTrial', ...
            'bpodTrialStart', 'zaberStepsUp', 'zaberStepsDwn', 'cameras', ...
            'zaberPerTrial', 'cameraPerTrial', ...
            'chan', 'eventMarkerDur' ...
            );
    end
    
    % For phy event_plugin_hh
    ks2Folder = dir(fullfile(fullfile(imecFolders(f).folder, imecFolders(f).name), 'imec*_ks2'));
    fullFileNamePhy = fullfile(ks2Folder.folder, ks2Folder.name, matName);
    copyfile(fullFileNameDJ, fullFileNamePhy);
    
    fprintf('%s saved!\n', fullFileNameDJ);
    fprintf('%s saved!\n', fullFileNamePhy);
    
end

function txtFile = getTxtFileName(sessionDir, chan, duration, chanTypeOverride)
if nargin < 4
    chanTypeOverride = '.XD';
end
% Parse non-integer duration
if mod(duration,1)
    durationStr = sprintf('%gp%g', floor(duration), mod(duration,1)*10);
else
    durationStr = sprintf('%g', duration);
end

txtFile = dir(fullfile(sessionDir, sprintf('*%s*%g_%s.adj.txt', chanTypeOverride, chan, durationStr)));  % Try TPrime adjusted first
if isempty(txtFile)  % If no adj.txt, try raw txt
    txtFile = dir(fullfile(sessionDir, sprintf('*%s*%g_%s.txt', chanTypeOverride, chan, durationStr)));
    fprintf('No _adj.txt found, using non-adjusted version');
end    
txtFile = fullfile(txtFile.folder, txtFile.name);
    
