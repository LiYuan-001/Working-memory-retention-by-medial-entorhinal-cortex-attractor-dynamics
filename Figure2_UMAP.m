
% plot figure for figure 2
% this code runs umap for the figure 8 maze data
% each time umap outcome in the 3-D coordinates might be different, but the
% relative distribution is similar
% Li YUAN, UCSD

clear all
close all

codePath = 'C:\Users\Li\Documents\MATLAB\DualProbe\CellProperty\Local\MEC_manuscript_code';
dataFolder = fullfile(codePath,'Data','1146_20240710');
addpath(genpath(codePath));

p.speed = 15; % unit: cm/s
p.timeBin = 200/10^3;

p.savePlot = 1;
p.writeToFile = 1;

% do umap for the run part maze and label with different way
% plot umap for this sequence
min_dist = 0.2;
spread = 2.0;
verbose = 'none';
method = 'mex';
metric = 'cosine';
see_training = false;
%  metric = 'euclidean';
n_neighbors = 20;
n_components = 3;
% -------------------------------------------------------------------------

% Read in input information
blockName = {'on10_1','off10_1','on30_1','off30_1','on10_2','off10_2','on30_2','off30_2'};
ana_block = [1:8];

ECLR_color = [0,0,1;1,0,1;1,0,0;0,1,1]; % blue, magenta, red, cyan

% syncTs = load(fullfile(sessInfo(i).NIDQ,'Processed', 'syncTime.mat'));
% load synchrinized pos t
load(fullfile(dataFolder, 'posT_Sync_NIDQ.mat'));
% load spike mat file
spkMatName = sprintf('%s%d%s','spkMat-',p.timeBin*10^3,'ms.mat');
spkMatFile = fullfile(dataFolder,spkMatName);
spkMat = load(spkMatFile);
% load cell int pyr label in MEC
cellTypeFile = fullfile(dataFolder,'mec_clusterType.mat');
load(cellTypeFile);

clusterNum = mec_clusterType.clusterNum;
pyrLabel = mec_clusterType.labelInd == 1;
clusterNum_Pyr = sum(pyrLabel);


pathData = loadPath(dataFolder,blockName);
% assign new position time to pos based on synced timestamps

for j = 1:length(blockName)

    spkTrain = spkMat.spkTrainMatrix.imec1.(blockName{j}).spkTrain_Raw;
    spkTrain_gauss = spkTrain;
    for k = 1:size(spkTrain,1)
        spkTrain_gauss(k,:) = gaussfilt(1:length(spkTrain(k,:)),spkTrain(k,:),2);
    end
    timeBinTemp = (spkMat.spkTrainMatrix.pos.(blockName{j}).timeBin(1:end-1)+spkMat.spkTrainMatrix.pos.(blockName{j}).timeBin(2:end))/2;

    posX = zeros(size(timeBinTemp));
    posY = zeros(size(timeBinTemp));
    posV = zeros(size(timeBinTemp));
    runLabel = zeros(size(timeBinTemp));
    delayTime = zeros(size(timeBinTemp));
    delayLabel = zeros(size(timeBinTemp));
    delayLabel_10s = zeros(size(timeBinTemp));
    delayLabel_block1 = zeros(size(timeBinTemp));
    delayLabel_block2 = zeros(size(timeBinTemp));
    delayLabel_block3 = zeros(size(timeBinTemp));
    delayLabel_on = zeros(size(timeBinTemp));
    delayLabel_off = zeros(size(timeBinTemp));


    delayFile = fullfile(dataFolder,'Position',blockName{j}, 'Fig8DelayZonePos.mat');
    load(delayFile);
    pdata = pathData.(blockName{j});

    pdata.t = posT_Sync_NIDQ.(blockName{j});

    % y2 = interp1(x1, y1, x2);
    posX = interp1(pdata.t',pdata.x',timeBinTemp);
    posY = interp1(pdata.t',pdata.y',timeBinTemp);
    posV = interp1(pdata.t',pdata.v,timeBinTemp);
    fastSpeedInd = posV >= p.speed;


    % take each trial to label delay period
    trialNum = length(Fig8DelayZonePos.delayPos1.startT_Sync);

    for n = 1:trialNum

        delayInd = timeBinTemp > (Fig8DelayZonePos.delayPos1.startT_Sync(n)+2) & timeBinTemp <= Fig8DelayZonePos.delayPos1.endT_Sync(n);
        % split the 30 s to 10 s intervals
        if contains(blockName{j},'30')
            block1_Ind = timeBinTemp > (Fig8DelayZonePos.delayPos1.startT_Sync(n)+2) & timeBinTemp <= (Fig8DelayZonePos.delayPos1.startT_Sync(n)+10);
            block2_Ind = (timeBinTemp > Fig8DelayZonePos.delayPos1.startT_Sync(n)+10) & timeBinTemp <= (Fig8DelayZonePos.delayPos1.startT_Sync(n)+20);
            block3_Ind = (timeBinTemp > Fig8DelayZonePos.delayPos1.startT_Sync(n)+20) & timeBinTemp <= (Fig8DelayZonePos.delayPos1.startT_Sync(n)+30);
            delayLabel_block1(block1_Ind) = 1;
            delayLabel_block2(block2_Ind) = 1;
            delayLabel_block3(block3_Ind) = 1;
        else
            block1_Ind = timeBinTemp > (Fig8DelayZonePos.delayPos1.startT_Sync(n)+2) & timeBinTemp <= (Fig8DelayZonePos.delayPos1.startT_Sync(n)+10);
            delayLabel_10s(block1_Ind) = 1;
        end
        delayTime(delayInd) = timeBinTemp(delayInd) - Fig8DelayZonePos.delayPos1.startT_Sync(n);
        delayLabel(delayInd) = 1;

        if contains(blockName{j},'on')
            delayLabel_on(delayInd) = 1;
        else
            delayLabel_off(delayInd) = 1;
        end
        runLabel(fastSpeedInd) = 1;
    end

    umap_DelayInterval_LR.(blockName{j}).posX = posX;
    umap_DelayInterval_LR.(blockName{j}).posY = posY;
    umap_DelayInterval_LR.(blockName{j}).posV = posV;

    umap_DelayInterval_LR.(blockName{j}).delayTime = delayTime;
    umap_DelayInterval_LR.(blockName{j}).delayLabel = delayLabel;
    umap_DelayInterval_LR.(blockName{j}).delayLabel_10s = delayLabel_10s;
    umap_DelayInterval_LR.(blockName{j}).delayLabel_block1 = delayLabel_block1;
    umap_DelayInterval_LR.(blockName{j}).delayLabel_block2 = delayLabel_block2;
    umap_DelayInterval_LR.(blockName{j}).delayLabel_block3 = delayLabel_block3;
    umap_DelayInterval_LR.(blockName{j}).delayLabel_on = delayLabel_on;
    umap_DelayInterval_LR.(blockName{j}).delayLabel_off = delayLabel_off;

    umap_DelayInterval_LR.(blockName{j}).runLabel = runLabel;
    umap_DelayInterval_LR.(blockName{j}).spkTrain_gauss = spkTrain_gauss;

end


% get the values required for plotting
%%
spkMat_run = [];
posX_allPos = [];
posY_allPos = [];
runLabel_all = [];
pos_run_Ind = [];
delayTimeAll = [];
delayLabel_Ind = [];
delayLabel_10s_Ind = [];
delayLabel_block1_Ind = [];
delayLabel_block2_Ind = [];
delayLabel_block3_Ind = [];
delayLabel_on_Ind = [];
delayLabel_off_Ind = [];

for j = 1:length(blockName)
    runLabel = umap_DelayInterval_LR.(blockName{j}).runLabel;
    runLabel_all = [runLabel_all,runLabel];
    pos_run_Ind = [pos_run_Ind,runLabel];


    delayTimeAll = [delayTimeAll,umap_DelayInterval_LR.(blockName{j}).delayTime];
    delayLabel_Ind = [delayLabel_Ind,umap_DelayInterval_LR.(blockName{j}).delayLabel];
    delayLabel_block1_Ind = [delayLabel_block1_Ind,umap_DelayInterval_LR.(blockName{j}).delayLabel_block1];
    delayLabel_block2_Ind = [delayLabel_block2_Ind,umap_DelayInterval_LR.(blockName{j}).delayLabel_block2];
    delayLabel_block3_Ind = [delayLabel_block3_Ind,umap_DelayInterval_LR.(blockName{j}).delayLabel_block3];

    delayLabel_10s_Ind = [delayLabel_10s_Ind,umap_DelayInterval_LR.(blockName{j}).delayLabel_10s];


    if contains(blockName{j},'on')
        delayLabel_on_Ind = [delayLabel_on_Ind,umap_DelayInterval_LR.(blockName{j}).delayLabel_on];
        delayLabel_off_Ind = [delayLabel_off_Ind,zeros(size(umap_DelayInterval_LR.(blockName{j}).delayLabel_on))];
    else
        delayLabel_off_Ind = [delayLabel_off_Ind,umap_DelayInterval_LR.(blockName{j}).delayLabel_off];
        delayLabel_on_Ind = [delayLabel_on_Ind,zeros(size(umap_DelayInterval_LR.(blockName{j}).delayLabel_off))];
    end


    posX_allPos = [posX_allPos,umap_DelayInterval_LR.(blockName{j}).posX];
    posY_allPos = [posY_allPos,umap_DelayInterval_LR.(blockName{j}).posY];
    spkMat_run = [spkMat_run,umap_DelayInterval_LR.(blockName{j}).spkTrain_gauss];
end

h = figure(1);
h.Position = [100,100,1600,600];
% plot maze pos in color
subplot(2,5,1)
colorx = (posX_allPos - min(posX_allPos)) / (max(posX_allPos) - min(posX_allPos));
colory = (posY_allPos - min(posY_allPos)) / (max(posY_allPos) - min(posY_allPos));
RGB = [colorx',colory',ones(size(colory'))*0.7];
scatter(posX_allPos, posY_allPos, 10, RGB, 'filled');
hold on
scatter(posX_allPos(delayLabel_on_Ind==1), posY_allPos(delayLabel_on_Ind==1), 30, [0.4,0.4,0.4], 'filled');
set(gca,'YDir','Reverse')
TITLE1 = 'Figure 8 maze';
TITLE2 = 'black delay on';
title({TITLE1;TITLE2},'Interpreter','None');

subplot(2,5,6)
colorx = (posX_allPos - min(posX_allPos)) / (max(posX_allPos) - min(posX_allPos));
colory = (posY_allPos - min(posY_allPos)) / (max(posY_allPos) - min(posY_allPos));
RGB = [colorx',colory',ones(size(colory'))*0.7];
scatter(posX_allPos, posY_allPos, 10, RGB, 'filled');
hold on
scatter(posX_allPos(delayLabel_off_Ind==1), posY_allPos(delayLabel_off_Ind==1), 30, [0,0,0], 'filled');
set(gca,'YDir','Reverse')
TITLE1 = 'Figure 8 maze';
TITLE2 = 'gray delay off';
title({TITLE1;TITLE2},'Interpreter','None');


% Calculate and plot umap for the whole Figure 8 maze
% note that each time running umap, the umap might rotate, th relative
% relation remain similar
% rotate mannually to get best angle to see the distribution
subplot(2,5,2)
% all maze
dataSet = spkMat_run(pyrLabel,:)';
dataSet_z = zscore(dataSet,0,1);
[reduction, umap, clusterIds, extras] = run_umap(dataSet_z,...
    'min_dist',min_dist,'metric',metric,'spread',spread,'n_components',n_components,'n_neighbors',n_neighbors, 'verbose', verbose);

reduction_run = reduction((delayLabel_on_Ind + delayLabel_off_Ind) ==0,:);
RGB_run = RGB((delayLabel_on_Ind+ delayLabel_off_Ind) ==0,:);
reduction_on = reduction(delayLabel_on_Ind==1,:);
reduction_off = reduction(delayLabel_off_Ind==1,:);
reduction_on_rgb = [0.5,0.5,0.5];
reduction_off_rgb = [0,0,0];

scatter3(reduction_run(:,1), reduction_run(:,2), reduction_run(:,3), 2, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
hold on
scatter3(reduction_on(:,1), reduction_on(:,2), reduction_on(:,3), 2, reduction_on_rgb, 'filled');
scatter3(reduction_off(:,1), reduction_off(:,2), reduction_off(:,3), 2, reduction_off_rgb, 'filled');
TITLE1 = 'all maze umap';
TITLE2 = sprintf('Manually roate to get best angle to view');
title({TITLE1;TITLE2},'Interpreter','None');

% rotate mannually to get best angle to see the distribution
subplot(2,5,7)
scatter3(reduction_run(:,1), reduction_run(:,2), reduction_run(:,3), 2, RGB_run, 'filled');
hold on
scatter3(reduction_on(:,1), reduction_on(:,2), reduction_on(:,3), 2, reduction_on_rgb, 'filled');
scatter3(reduction_off(:,1), reduction_off(:,2), reduction_off(:,3), 2, reduction_off_rgb, 'filled');

TITLE1 = 'all maze umap rotated';
TITLE2 = sprintf('%s%1.2f%s%1.1f%s%s%s%d','min_dist:',min_dist,' spread:',spread,' metric:',metric,' n_neighbors:',n_neighbors);
title({TITLE1;TITLE2},'Interpreter','None');

[az0, el0] = view;   % save the original 3D view
view(az0+180, el0)   % rotate 180 around z, keep 3D angle

% treadmill on
% rotate mannually to get best angle to see the distribution
subplot(2,5,3)
valIdx = find((delayLabel_10s_Ind & delayLabel_on_Ind)==1);
colorVal = delayTimeAll(valIdx);
scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 2, colorVal, 'filled'); % 50 = marker size, 'filled' = solid circles
colormap(jet)
% colorbar
clim([0 30])
TITLE1 = 'delay on 10 s';
TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
title({TITLE1;TITLE2},'Interpreter','None');

% whole 30 s
% rotate mannually to get best angle to see the distribution
subplot(2,5,4)
delay_30_Ind = (delayLabel_block1_Ind + delayLabel_block2_Ind + delayLabel_block3_Ind) > 0;
valIdx = find((delay_30_Ind & delayLabel_on_Ind)==1);
colorVal = delayTimeAll(valIdx);
scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 2, colorVal, 'filled'); % 50 = marker size, 'filled' = solid circles
colormap(jet)
% colorbar
clim([0 30])
TITLE1 = 'delay on 0 - 30 s';
TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
title({TITLE1;TITLE2},'Interpreter','None');
[az0, el0] = view;   % save the original 3D view
% view(az0, el0-30)  

% whole 30 s
% rotate mannually to get best angle to see the distribution
subplot(2,5,5)
delay_30_Ind = (delayLabel_block1_Ind + delayLabel_block2_Ind + delayLabel_block3_Ind) > 0;
valIdx = find((delay_30_Ind & delayLabel_on_Ind)==1);
colorVal = delayTimeAll(valIdx);
scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 2, colorVal, 'filled'); % 50 = marker size, 'filled' = solid circles
colormap(jet)
% colorbar
clim([0 30])
TITLE1 = 'different projection';
TITLE2 = '';
title({TITLE1;TITLE2},'Interpreter','None');

[az0, el0] = view;   % save the original 3D view
view(az0+180, el0)   

% treadmill off
% rotate mannually to get best angle to see the distribution
subplot(2,5,8)
valIdx = find((delayLabel_10s_Ind & delayLabel_off_Ind)==1);
colorVal = delayTimeAll(valIdx);
scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 2, colorVal, 'filled');
colormap(jet)
% colorbar
clim([0 30])
TITLE1 = 'delay off 10 s';
TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
title({TITLE1;TITLE2},'Interpreter','None');


% whole 30 s off
% rotate mannually to get best angle to see the distribution
subplot(2,5,9)
delay_30_Ind = (delayLabel_block1_Ind + delayLabel_block2_Ind + delayLabel_block3_Ind) > 0;
valIdx = find((delay_30_Ind & delayLabel_off_Ind)==1);
RGB_run = zeros(length(valIdx),3);
colorVal = delayTimeAll(valIdx);
scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 2, colorVal, 'filled'); % 50 = marker size, 'filled' = solid circles
colormap(jet)
% colorbar
clim([0 30])
TITLE1 = 'delay off 0 - 30 s';
TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
title({TITLE1;TITLE2},'Interpreter','None');

% plot 10-s delay in L vs R return
subplot(2,5,10)
clim([0 30])
colorbar

fprintf('Finished figure 2\n');

