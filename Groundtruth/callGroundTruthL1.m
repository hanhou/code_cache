%% == Mannual selection ==
imecDir = uigetdir('F:\catGT\HH102');
groundTruthL1(imecDir);

%% == Batch processing ==
rootDirs = {'F:\catGT\HH102\', 'E:\catGT'};
includeList = {'imec1'};

for dd = 1:length(rootDirs)
    thisDir = rootDirs{dd};
    allMeta = dir(fullfile(thisDir, '**', '*.ap.meta'));

    for fN = 1:length(allMeta)
        imecDir = allMeta(fN).folder;
        if ~isempty(includeList) && ~contains(imecDir, includeList)
            fprintf('\nSkipped: %s\n', imecDir);
            continue
        end
        
        fprintf('\n===== Batch processing =====\n%s\n', imecDir);
        try
            groundTruthL1(imecDir);
        catch
            disp('ERROR!!')
        end
        close all;
    end
end


%% 
figs = findall(0, 'type', 'figure');
for f=1:length(figs)
    figure(f); xlim([-2 10])
end
