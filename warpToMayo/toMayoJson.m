%% Based on Dave's code
%% get xyz_picks.json
[fN, dirN] = uigetfile('F:\Z1\*.csv'); 
% pts = readBigWarpLandmarks(fullfile(dirN, fN));
pts = xlsread(fullfile(dirN,fN));
pts = pts(:, 3:end);

% lastPt=find(pts(:,1)==Inf,1,'first');
lastPt=find(isnan(pts(:,1)),1,'first');
pts=flip(round(pts(lastPt:end,4:6)));
xyz_pick=[pts(:,1)-5739 5400-pts(:,3) 332-pts(:,2)];
csStruct=struct('xyz_picks',xyz_pick);
fid = fopen([dirN 'xyz_picks.json'], 'w');
fwrite(fid, jsonencode(csStruct), '*char');
fclose(fid);
disp('Json saved.')