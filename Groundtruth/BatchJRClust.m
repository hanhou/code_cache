rootDirs = {'F:\catGT\HH102\', 'E:\catGT'};
includeList = {};

%% Batch bootstrap and sorting
for dd = 1:length(rootDirs)
    thisDir = rootDirs{dd};
    allMeta = dir(fullfile(thisDir, '**', '*.ap.meta'));
    
    for fN = 1:length(allMeta)
        metaName = allMeta(fN).name;
        metaFolder = allMeta(fN).folder;
        if ~isempty(includeList) && ~contains(metaName, includeList)
            fprintf('\nSkipped: %s\n', metaName);
            continue
        end
        fprintf('\n======= Batch JRClust =======\n%s\n', metaName);
        
        try
            % -- Bootstrap --
            eval(sprintf('jrc bootstrap %s -noconfirm -advanced -noedit', fullfile(metaFolder, metaName)));
            prmFile = fullfile(metaFolder, [metaName(1:end-5) '.prm']);
            
            % -- Detect and sorting --
            eval(sprintf('jrc detect-sort %s', prmFile));
        catch
            fprintf('ERROR!!!')
        end
        
        
    end
end
