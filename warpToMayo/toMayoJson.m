%% Based on Dave's code
%% Convert .csv from BigWarp to xyz_picks.json for Mayo's GUI

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
    name = regexp(fN,'(.*imec\d)(_shank)?(\d)?.csv', 'tokens');
    if length(name{1}) > 1
        newName = sprintf('%s_xyz_picks%s%g.json', name{1}{1}, name{1}{2}, str2num(name{1}{3})+1);
    else
        newName = sprintf('%s_xyz_picks.json', name{1}{1});
    end
    
    % Save file
    saveFileName = [dirN newName];
    fid = fopen(saveFileName, 'w');
    fwrite(fid, jsonencode(csStruct), '*char');
    fclose(fid);
    disp(['Json saved: ' saveFileName])
    
    % To-do: copy to corresponding alf folders

end

