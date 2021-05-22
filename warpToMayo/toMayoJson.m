%% Based on Dave's code
%% get xyz_picks.json
[fNs, dirN] = uigetfile('E:\MAPAlignment\*.*', 'Multiselect', 'on'); 
% pts = readBigWarpLandmarks(fullfile(dirN, fN));
if ~iscell(fNs) fNs = {fNs} ;end;

for ff = 1:length(fNs)
    fN = fNs{ff};
    pts = xlsread(fullfile(dirN,fN));
    pts = pts(:, 3:end);

    % lastPt=find(pts(:,1)==Inf,1,'first');
    lastPt=find(isnan(pts(:,1)),1,'first');
    pts=flip(round(pts(lastPt:end,4:6)));
    xyz_pick=[pts(:,1)-5739 5400-pts(:,3) 332-pts(:,2)];
    csStruct=struct('xyz_picks',xyz_pick);
    
    saveFileName = [dirN sprintf('%s_xyz_picks.json', fN(1:end-4))];
    fid = fopen(saveFileName, 'w');
    fwrite(fid, jsonencode(csStruct), '*char');
    fclose(fid);
    disp(['Json saved: ' saveFileName])
end