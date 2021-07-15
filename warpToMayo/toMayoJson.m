%% Based on Dave's code
%% Convert .csv from BigWarp to xyz_picks.json for Mayo's GUI
%% Han Hou 2021

[fNs, dirN] = uigetfile('f:\Z1\*.*', 'Multiselect', 'on'); 
if ~iscell(fNs) fNs = {fNs}; end;

for ff = 1:length(fNs)
    % Get data
    fN = fNs{ff};
    pts = xlsread(fullfile(dirN,fN));
    pts = pts(:, end-5:end);
    
    % Convert CCF coordinate to Mayo's coordinate
    lastPt=find(isnan(pts(:,1)),1,'first');
    pts=flip(round(pts(lastPt:end,4:6)));
    xyz_pick=[pts(:,1)-5739 5400-pts(:,3) 332-pts(:,2)];
    csStruct=struct('xyz_picks',xyz_pick);
    
    % Convert file name to Mayo's convention (..._imecX_xyz_picks_shankX.json)
    name = regexp(fN,'(.*\d{4})_(\d)_?(\d)?', 'tokens');  
    if isempty(name)
        sprintf('Wrong file name: %s!!\n', fN);
        continue
    elseif isempty(name{1}{3})
        imec = str2double(name{1}{2})-1;   % Mayo uses imec0, 1
        newName = sprintf('%s_imec%g_xyz_picks.json', name{1}{1}, imec);  
    else
        imec = str2double(name{1}{2})-1;   % Mayo uses imec0, 1
        shank = str2double(name{1}{3});     % Mayo uses shank 1,2,3,4
        newName = sprintf('%s_imec%g_xyz_picks_shank%g.json', name{1}{1}, imec, shank); 
    end
    
    % Save file
    saveFileName = [dirN newName];
    fid = fopen(saveFileName, 'w');
    fwrite(fid, jsonencode(csStruct), '*char');
    fclose(fid);
    disp(['Json saved: ' saveFileName])
    
    % To-do: copy to corresponding alf folders

end

