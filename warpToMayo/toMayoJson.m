%% Based on Dave's code
%% Convert .csv from BigWarp to xyz_picks.json for Mayo's GUI
%% Han Hou 2021

% Name format
% NP1.0 and 2.1: "landmarks_HH08_20210819_1.csv"
% NP2.4:         "landmarks_HH08_20210824_1_1.csv"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fNs, dirN] = uigetfile('f:\*.*', 'Multiselect', 'on'); 
alf_root_folder = 'i:\catGT\HH08\';  % If not empty, automatically copy the .json file to corresponding alf folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(fNs) fNs = {fNs}; end;

alf_list = dir(fullfile(alf_root_folder, '**\alf'));  %get list of files and folders in any subfolder
alf_list = {alf_list([alf_list(:).isdir]).folder}';

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
    name = regexp(fN,'(.*\d{4})_(\d)_?(\d|npx.*)?', 'tokens');    % Account for possible "npx" suffix (only for TrackFinder)
    if isempty(name)
        fprintf('Wrong file name: %s!!\n', fN);
        continue
    elseif isempty(name{1}{3}) || contains(name{1}{3}, 'npx')  % NP1.0 or NP2.1
        imec = str2double(name{1}{2})-1;   % Mayo uses imec0, 1
        newName = sprintf('%s_imec%g_xyz_picks.json', name{1}{1}, imec);  
    else  % NP2.4
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
    
    % Copy to corresponding alf folders
    if ~isempty(alf_list)
        % Locate alf folder
        alf_pattern = sprintf('%s_g\\d_imec%g', name{1}{1}(end-7:end), imec);  % e.g., '20210416_g\d_imec0'
        for aa = 1:length(alf_list)
            alf = alf_list{aa};
            if ~isempty(regexp(alf, alf_pattern, 'Once'))
                % Copy file
                copyfile(saveFileName, alf);
                fprintf('  Copyed to %s!\n', alf);
                break
            end
        end
    end
end

