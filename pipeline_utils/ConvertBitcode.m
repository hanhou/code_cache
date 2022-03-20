% Based on Dave's code
function ConvertBitcode(sessionDir)
%% Settings
if nargin == 0
    sessionDir=uigetdir('I:\catGT\');
end

% In line with my pybpod code
bitCodeDigits = 20;

chan.protocol = 2;  % Protocol dig marker channel
eventMarkerDur.bitcodeEachbit = '10';
eventMarkerDur.bitcodeFirst = '20';  % The multiplier of the first bit code width
eventMarkerDur.bitcodeFirstMultiplier = 2;

eventMarkerDur.goCue = '1';
eventMarkerDur.reward = '30';
eventMarkerDur.ITI = '40';

chan.behavior = 2;  % Behavior dig marker channel (I moved all behavior channel to bitcode channel in Apr 21 version)
eventMarkerDur.choiceL = '2';  
eventMarkerDur.choiceR = '3';

chan.bpodstart = 1;
chan.zaber = 4;
chan.cameras = [5, 6, 7];
chan.cameraNameInDJ = {'Camera 1', 'Camera 0', 'Camera 2'};   % (From bottom to up): Chameleon3 CM3-U3-13Y3M-CS (FLIR); 
                                                              % 300 Hz Bottom face, 300 Hz Side face, 100 Hz Body
                                                              
chan.leftLaser = 1; % XA
chan.rightLaser = 2;  % XA
chan.leftLick = 3; % XA
chan.rightLick = 4;  % XA

%% get Ephys Bitcode
% Start of a trial (onset of my bitcode is indicated by twice of the bitcode width)
sTrig = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.bitcodeFirst)); 
iti = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.ITI));

bitsAll = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.bitcodeEachbit));
goCue = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.goCue));
reward = dlmread(getTxtFileName(sessionDir, chan.protocol, eventMarkerDur.reward));
choiceL = dlmread(getTxtFileName(sessionDir, chan.behavior, eventMarkerDur.choiceL));
choiceR = dlmread(getTxtFileName(sessionDir, chan.behavior, eventMarkerDur.choiceR));
choiceR = setdiff(choiceR, choiceL);   % 2.5 * (1-20%) = 2.0!!! choice R also includes all choice L...

bpodTrialStart = dlmread(getTxtFileName(sessionDir, chan.bpodstart, '0'));
zaberStepsUp = dlmread(getTxtFileName(sessionDir, chan.zaber, '0'));
zaberStepsDwn = dlmread(getTxtFileName(sessionDir, chan.zaber, '0', '.iXD'));
zaberStepsAll = sort([zaberStepsUp; zaberStepsDwn]);

for i = 1:length(chan.cameras)
    cameras{i} = dlmread(getTxtFileName(sessionDir, chan.cameras(i), '0'));
end

bitcode = zeros(length(iti), bitCodeDigits);   % Use length iti to make sure (the last) trial has ended.

% All licks
lickLAll = dlmread(getTxtFileName(sessionDir, chan.leftLick, '0', 'XA'));
lickRAll = dlmread(getTxtFileName(sessionDir, chan.rightLick, '0', 'XA'));

% Laser pulses
try
    laserLAll = dlmread(getTxtFileName(sessionDir, chan.leftLaser, '0', 'XA'));
catch % In case there are no laser pulses
    laserLAll = [];
end

try
    laserRAll = dlmread(getTxtFileName(sessionDir, chan.rightLaser, '0', 'XA'));
catch % In case there are no laser pulses
    laserRAll = [];
end

% headings in datajoint pipeline, ephys.TrialEventType
headings = {'bitcodestart', 'go', 'choice', 'choice', 'reward', 'trialend', 'bpodstart', 'zaberready'};
[STRIG_, GOCUE_, CHOICEL_, CHOICER_, REWARD_, ITI_, BPOD_START_, ZABER_IN_POS_] = deal(1,2,3,4,5,6,7,8);

% Use iti as trial marker to exclude truncated trials
digMarkerPerTrial = nan(length(iti), 8);  % [STrig, goCue, choiceL, choiceR, reward, ITI]
zaberPerTrial = {};
cameraPerTrial = {};
lickLPerTrial = {};
lickRPerTrial = {};

% digMarkerPerTrial(:,[STRIG_, GOCUE_, ITI_]) = [sTrig, goCue, iti];  % Must exists 
digMarkerPerTrial(:,ITI_) = iti;  % Must exists 
countFakeZaberStep = 0;

for i = 1:length(iti)  % Fill in digMarkerPerTrial
    thisBitCodeStart = sTrig(find(sTrig < iti(i), 1, 'last'));  % Deal with truncated trials (only search backward from each iti)
    if i < length(iti)
        nextBitCodeStart = sTrig(find(sTrig < iti(i + 1), 1, 'last'));
    else
        nextBitCodeStart = Inf;
    end

    digMarkerPerTrial(i, BPOD_START_) = bpodTrialStart(find(bpodTrialStart < iti(i), 1, 'last'));   % The closest bpodTrialStart before this ITI
    digMarkerPerTrial(i, STRIG_) = thisBitCodeStart;
    digMarkerPerTrial(i, GOCUE_) = goCue(thisBitCodeStart < goCue & goCue < iti(i));
        
    % Parse bitcode
    bitHighThis = bitsAll(thisBitCodeStart < bitsAll & bitsAll < iti(i));
    bitHighPositionThis = round((bitHighThis - thisBitCodeStart - ...
        (eventMarkerDur.bitcodeFirstMultiplier - 1) * str2double(eventMarkerDur.bitcodeEachbit) * 1e-3) ...
                                 / (2 * str2double(eventMarkerDur.bitcodeEachbit) * 1e-3));
    
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
    if isempty(zaberPerTrial{i})
        % For some rare cases where zaber feedback were missing, fake the
        % last zaber step
        fakeZaberStep = thisBitCodeStart + (2 * str2double(eventMarkerDur.bitcodeEachbit) * bitCodeDigits ...
            + str2double(eventMarkerDur.bitcodeFirst))/1000 + 0.136; 
        zaberPerTrial{i} = fakeZaberStep;
        digMarkerPerTrial(i, ZABER_IN_POS_) = fakeZaberStep;
        countFakeZaberStep = countFakeZaberStep + 1;
    else
        digMarkerPerTrial(i, ZABER_IN_POS_) = max(zaberPerTrial{i});   % The last zaber pulse of this trial
    end
    
    % all licks (from this trial start to the next trial start, including iti after this trial; just in case where we don't retract lickports)
    lickLPerTrial{i} = lickLAll(thisBitCodeStart < lickLAll & lickLAll < nextBitCodeStart);
    lickRPerTrial{i} = lickRAll(thisBitCodeStart < lickRAll & lickRAll < nextBitCodeStart);
    
    % photo stimulation (from this trial start to next trial start, including iti after this trial; Note: for laserXPerTrial{i}, we expect to see behavioral effects on trial i+1)
    laserLPerTrial{i} = laserLAll(thisBitCodeStart < laserLAll & laserLAll < nextBitCodeStart);
    laserRPerTrial{i} = laserRAll(thisBitCodeStart < laserRAll & laserRAll < nextBitCodeStart);
     
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

if countFakeZaberStep > 0
    fprintf('Zaber steps faked for %g trials!\n', countFakeZaberStep)
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
            'bitcode', 'bitCodeS', 'goCue', 'sTrig', 'reward', 'choiceL', 'choiceR', 'choiceAll', 'iti', ... 
            'digMarkerPerTrial', 'headings', ...
            'bpodTrialStart', 'zaberStepsUp', 'zaberStepsDwn', 'cameras', ...
            'chan', 'eventMarkerDur', ...
            'zaberPerTrial', 'cameraPerTrial', ...
            'lickLAll', 'lickRAll', 'lickLPerTrial', 'lickRPerTrial', ...
            'laserLPerTrial', 'laserRPerTrial', ...
            'trialNum');
        fprintf('Trial Number fixed!!\n')
    else
        save(fullFileNameDJ, ...
            'bitcode', 'bitCodeS', 'goCue', 'sTrig', 'reward', 'choiceL', 'choiceR', 'choiceAll', 'iti', ...
            'digMarkerPerTrial', 'headings',...
            'bpodTrialStart', 'zaberStepsUp', 'zaberStepsDwn', 'cameras', ...
            'chan', 'eventMarkerDur', ...
            'zaberPerTrial', 'cameraPerTrial', ...
            'lickLAll', 'lickRAll', 'lickLPerTrial', 'lickRPerTrial', ...
            'laserLPerTrial', 'laserRPerTrial' ...
            );
    end
    
    % For phy event_plugin_hh
    ks2Folder = dir(fullfile(fullfile(imecFolders(f).folder, imecFolders(f).name), 'imec*_ks2'));
    fullFileNamePhy = fullfile(ks2Folder.folder, ks2Folder.name, matName);
    copyfile(fullFileNameDJ, fullFileNamePhy);
    
    fprintf('%s saved!\n', fullFileNameDJ);
    fprintf('%s saved!\n', fullFileNamePhy);
    
end

function output = getTxtFileName(sessionDir, chan, durationStr, chanTypeOverride)
if nargin < 4
    chanTypeOverride = '.XD';
end

txtFile = dir(fullfile(sessionDir, sprintf('*%s*%g_%s.adj.txt', chanTypeOverride, chan, durationStr)));  % Try TPrime adjusted first
if isempty(txtFile)  % If no adj.txt, try raw txt
    txtFile = dir(fullfile(sessionDir, sprintf('*%s*%g_%s.txt', chanTypeOverride, chan, durationStr)));
    fprintf('No _adj.txt found, using non-adjusted version\n');
end

output = strings(1,length(txtFile));
for f = 1:length(txtFile)
    output(f) = fullfile(txtFile(f).folder, txtFile(f).name);
end
    
