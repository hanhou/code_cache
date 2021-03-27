% Based on Dave's code
clear all; close all; clc
%%
filepath='H:\data\MAP\dl62\dl62_012819_g0\';
%% get Ephys Bitcode
digital=dlmread([filepath filepath(end-15:end-1) '_imec1\' filepath(13:16) '_g0_tcat.imec.SY_384_0_0.txt'])
sTrig=digital;
digital=dlmread([filepath filepath(end-15:end-1) '_imec1\' filepath(13:16) '_g0_tcat.imec.SY_384_1_50.txt']);
eTrig=[sTrig(2:end); sTrig(end)+20];
goCue=digital;
digital=dlmread([filepath filepath(end-15:end-1) '_imec1\' filepath(13:16) '_g0_tcat.imec.SY_384_1_10.txt']);
bitcode=zeros(length(sTrig),20);
​
for i = 1:length(sTrig)
    bitHigh=digital(digital(:,1)>sTrig(i) & digital(:,1)<eTrig(i),:);
    
    preBitcode=cumsum(round(diff(bitHigh)/0.02));
    preBitcode=preBitcode(preBitcode<21);
    bitcode(i,preBitcode)=1;
    if i == 1
        bitCodeS=sprintf('%d', bitcode(i,:));
    else
        bitCodeS=[bitCodeS; sprintf('%d', bitcode(i,:))];
    end
end
​
%% get bitcode
clearvars -except bitcode bitCodeS goCue sTrig
%%
save([filepath filepath(end-15:end-1) '_imec1\bitcode.mat'])