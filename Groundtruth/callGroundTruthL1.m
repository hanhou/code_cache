%% == Mannual selection ==
% imecDir = uigetdir('F:\catGT\HH102');
% groundTruthL1(imecDir);

%% == Batch processing ==
rootDirs = {'F:\catGT\HH102\', 'E:\catGT'};
includeList = {};

for dd = 1:length(rootDirs)
    thisDir = rootDirs{dd};
    allMeta = dir(fullfile(thisDir, '**', '*.ap.meta'));

    try
        for fN = 1:length(allMeta)
            imecDir = allMeta(fN).folder;
            groundTruthL1(imecDir);
        end
    catch
        
    end
end
