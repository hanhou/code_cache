%% Batch convert txt files from catGT to bitcode.mat files for phy and ingestion
rootPath = 'F:\catGT';
allSubFolder = dir(fullfile(rootPath, '**\*.*'));
allSubFolder = allSubFolder([allSubFolder.isdir]);
allSubFolder = allSubFolder(arrayfun(@(A)contains(A.name,'catgt'), allSubFolder));

%%
for f = 1:length(allSubFolder)
    thisFolder = fullfile(allSubFolder(f).folder, allSubFolder(f).name);
    try
        ConvertBitcode(thisFolder);
    catch
        fprintf('!!! Errors in processing %s\n', thisFolder);
    end
end