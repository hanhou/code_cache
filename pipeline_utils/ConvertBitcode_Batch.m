%% Batch convert txt files from catGT to bitcode.mat files for phy and ingestion
rootPath = 'K:\';
allSubFolder = dir(fullfile(rootPath, '**\*.*'));
allSubFolder = allSubFolder([allSubFolder.isdir]);
allSubFolder = allSubFolder(arrayfun(@(A)contains(A.name,'catgt'), allSubFolder));
exclude = {'surface'};

%%
for f = 1:length(allSubFolder)
    thisFolder = fullfile(allSubFolder(f).folder, allSubFolder(f).name);
    if contains(thisFolder, exclude)
        continue
    end
    
    try
        ConvertBitcode(thisFolder);
    catch exception
        getReport(exception)
        fprintf('!!! Errors in processing %s\n', thisFolder);
    end
end