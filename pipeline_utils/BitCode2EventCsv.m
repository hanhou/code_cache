%% Batch convert GoCue.mat or BitCode.mat to Event.csv for spike curation in phy
rootPath = 'X:\Curation';
pattern = {'bitCode', 'GoCue'};
allMatFile = dir(fullfile(rootPath, '**\*.mat'));

%%
for f = 1:length(allMatFile)
    this = allMatFile(f);
    if ~contains(this.name, pattern, 'IgnoreCase', true)
        continue
    end
    
    thisMat = fullfile(this.folder, this.name);
    fprintf('Converting: %s..\n', thisMat);
    result = load(thisMat);
    try
        goCue = result.goCue;
    catch
        fprintf('  !! Data invalid !!')
        continue
    end
    
    writematrix(goCue, fullfile(this.folder, 'events.csv'));
    ksFolder = fullfile(this.folder, dir(fullfile(this.folder, '*_ks*')).name);
    copyfile(fullfile(this.folder, 'events.csv'), ksFolder);
    
 end