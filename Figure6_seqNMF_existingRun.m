% plot figure for figure 6
% this code plots exisiting sequence which is same as the figure 6
% there is another code calculate seqNMF from cell matrix. Each run of
% seqNMF generate slightly different results due to randomization, but the
% general outcome is similar
% Li YUAN, UCSD

close all
clear all
% % Turn off figure pop-ups only for this code block
% set(0, 'DefaultFigureVisible', 'on');

codePath = 'C:\Users\Li\Documents\MATLAB\DualProbe\CellProperty\Local\MEC_manuscript_code'; % change the folder to your local path
dataFolder = fullfile(codePath,'Data','1146_20240711');
addpath(genpath(codePath));

p.timeBin = 200/10^3;
p.L_time = 5; % unit in second
p.seqThresPct = 99;
p.consBinThres = 6;
p.w_min = 1;
p.min_pyrNum = 10; % here make it low to include more sequences, will increase to 15 for quant
p.pValThres = 0.05;
p.minCorrCell = 9; % because matalb corr doesnt calculate p if < 9

% linear map area
p.linROI = [33, 46];
p.linChoice = [55, 64];
p.xTickLin = [16 40 50 60 68];
p.xTickLabel = {'Return','Dl','St','Ch','Rw'};

p.distanceBin_return = [50,80,110,140];
p.distanceBin_stem = [50,80];
p.distIncre = 10;

p.savePlot = 1;
p.writeToFile = 0;

% Read in input information
blockName = {'on10_1','off10_1','on30_1','off30_1','on10_2','off10_2','on30_2','off30_2'};

map_All_imec1  = [];
map_On_imec1 = [];
map_Off_imec1  = [];

% load synchrinized pos t. 
% synchronized neuralynx to neuropixels time
load(fullfile(dataFolder,'posT_Sync_NIDQ.mat'));
% load cell putative principle label in MEC
cellTypeFile = fullfile(dataFolder, 'mec_clusterType.mat');
load(cellTypeFile);
% load single pike file
imec1_file = fullfile(dataFolder, 'pixelCluster_imec1.mat');
load(imec1_file);
% load speed
delaySpeedFile = fullfile(dataFolder,'Delay_Speed.mat');
load(delaySpeedFile);

% load sequences
fileName = sprintf('%s%d%s%d%s','Fig8DelaySeq-',ceil(p.L_time*10^3),'ms-',p.timeBin*10^3,'ms.mat');
seq_file = fullfile(dataFolder, fileName);
load(seq_file);

% load shuffled sequences
fileName = sprintf('%s%d%s%d%s','Fig8CenterSeq-Shuffle-',ceil(p.L_time*10^3),'ms-',p.timeBin*10^3,'ms.mat');
seq_file = fullfile(dataFolder, fileName);
load(seq_file);

% load extracted seq
fileName = sprintf('%s%d%s%d%s','seqNMF_DelayDetect_Event_Extract_peak-',ceil(p.L_time*10^3),'ms-',p.timeBin*10^3,'ms.mat');
seq_file = fullfile(dataFolder, fileName);
load(seq_file);

% load delayFire map
delayFile = fullfile(dataFolder,'Fig8DelayTimeMap_2Session.mat');
load(delayFile);

% get valid cell ind (putative principle cells)
prcLabel_imec1 = Fig8_seq_Delay_detect.imec1.pyrInd;
% avgRate_imec1 = mec_clusterType.avgRate;
depth_imec1 = mec_clusterType.clusterDepth;
cellLayer_imec1 = mec_clusterType.layerInd;
depthLayer_imec1 = mec_clusterType.posInd;
clusterID_imec1 = mec_clusterType.clusterID;
clusterNum = mec_clusterType.clusterNum;


% get average lin map for each block
spikeMap = cell(clusterNum,1);
spikeMap_On = cell(clusterNum,1);
spikeMap_Off = cell(clusterNum,1);
occupMap = cell(clusterNum,1);
occupMap_On = cell(clusterNum,1);
occupMap_Off = cell(clusterNum,1);

for j = 1:length(blockName)
    % load analyzed map per trial
    load(fullfile(dataFolder,'Position',blockName{j},'ratesByECLR_imec1.mat'));
    trialNum = length(ratesByECLR_imec1.ECLR);

    for n = 1:trialNum
        for k = 1:clusterNum
            spikeMap{k} = [spikeMap{k};ratesByECLR_imec1.spikemap{k,n}];
            occupMap{k} = [occupMap{k};ratesByECLR_imec1.occupmap{k,n}];

            if contains(blockName{j},'on')
                spikeMap_On{k} = [spikeMap_On{k};ratesByECLR_imec1.spikemap{k,n}];
                occupMap_On{k} = [occupMap_On{k};ratesByECLR_imec1.occupmap{k,n}];
            elseif contains(blockName{j},'off')
                spikeMap_Off{k} = [spikeMap_Off{k};ratesByECLR_imec1.spikemap{k,n}];
                occupMap_Off{k} = [occupMap_Off{k};ratesByECLR_imec1.occupmap{k,n}];
            else
                error('Block type wrong')
            end

        end
    end

end

for k = 1:clusterNum
    spikeMap_AllTrial = sum(spikeMap{k},1);
    occupMap_AllTrial = sum(occupMap{k},1);
    map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
    map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
    map_All_imec1 = [map_All_imec1;map_Temp_2];

    % % treadmill on
    % spikeMap_AllTrial = sum(spikeMap_On{k},1);
    % occupMap_AllTrial = sum(occupMap_On{k},1);
    % map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
    % map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
    % map_On_imec1 = [map_On_imec1;map_Temp_2];
    % 
    % % treadmill off
    % spikeMap_AllTrial = sum(spikeMap_Off{k},1);
    % occupMap_AllTrial = sum(occupMap_Off{k},1);
    % map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
    % map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
    % map_Off_imec1 = [map_Off_imec1;map_Temp_2];

end

% get each phase names (no delay etc)
blockName2 = {'on10','off10','on30','off30'};
for j = 1:length(blockName2)

    timeMap_Def1.(blockName2{j}) = [];
    timeMap_Def2.(blockName2{j}) = [];
    pre_Map_Def1.(blockName2{j}) = [];
    pre_Map_Def2.(blockName2{j}) = [];
    post_Map_Def1.(blockName2{j}) = [];
    post_Map_Def2.(blockName2{j}) = [];

    for k = 1:clusterNum
        timeMap_Def1.(blockName2{j}) = [timeMap_Def1.(blockName2{j});Fig8DelayTimeMap_2Session.imec1.(blockName2{j}).spikeRate1_Combined_Smooth{k}];
        timeMap_Def2.(blockName2{j}) = [timeMap_Def2.(blockName2{j});Fig8DelayTimeMap_2Session.imec1.(blockName2{j}).spikeRate2_Combined_Smooth{k}];

        pre_Map_Def1.(blockName2{j}) = [pre_Map_Def1.(blockName2{j});Fig8DelayTimeMap_2Session.imec1.(blockName2{j}).pre_spikeRate1_Combined_Smooth{k}];
        pre_Map_Def2.(blockName2{j}) = [pre_Map_Def2.(blockName2{j});Fig8DelayTimeMap_2Session.imec1.(blockName2{j}).pre_spikeRate2_Combined_Smooth{k}];

        post_Map_Def1.(blockName2{j}) = [post_Map_Def1.(blockName2{j});Fig8DelayTimeMap_2Session.imec1.(blockName2{j}).post_mapMean_Def1_2Session{k}];
        post_Map_Def2.(blockName2{j}) = [post_Map_Def2.(blockName2{j});Fig8DelayTimeMap_2Session.imec1.(blockName2{j}).post_mapMean_Def2_2Session{k}];
    end
end

% get weight dist in seqNMF for all pyr cells
shuffle_pyr = Fig8_seq_Center_detect_Shuffle.imec1.pyrInd;
if any(prcLabel_imec1 ~= shuffle_pyr)
    error('Error in cell IDs')
end

seqNum_Shuffle = Fig8_seq_Center_detect_Shuffle.K;
shuffleTimes = Fig8_seq_Center_detect_Shuffle.shuffleTimes;
shuffleWeight = zeros(sum(prcLabel_imec1),shuffleTimes*seqNum_Shuffle);
for n = 1:shuffleTimes
    w_shuffle = Fig8_seq_Center_detect_Shuffle.imec1.zpyr.W{n};
    %         w_shuffle(:,:,1:2)=0;
    %         w_shuffle(:,:,end-1:end)=0;
    w_shuffle2 = max(w_shuffle,[],3);
    shuffleWeight(:,(n-1)*seqNum_Shuffle+(1:seqNum_Shuffle)) = w_shuffle2;
end
% Sort each row in ascending order
sortedWeights = sort(shuffleWeight, 2);
% Determine the exact index for the desired percentile
N = size(shuffleWeight, 2);  % number of columns
NN = ceil(p.seqThresPct / 100 * N);  % index of percentile (round up)
% Extract the exact values
w_thres = sortedWeights(:, NN);
w_thres2 = w_thres;
w_thres2(w_thres ==0) = NaN;


% start taking cells from sequences
seqNum = Fig8_seq_Delay_detect.K;
seqWght_prc = Fig8_seq_Delay_detect.imec1.zpyr.W;
seqStrength = Fig8_seq_Delay_detect.imec1.zpyr.H_alltime;
binTime = Fig8_seq_Delay_detect.timeBin_all;
indSort_pyr = Fig8_seq_Delay_detect.imec1.zpyr.indSort;
% sort for all sequence visualization
seqWght_pyr_sort = seqWght_prc(indSort_pyr,:,:);


for k = 3
    seqName = sprintf('seq%d',k);

    % sort by this sequences for visualization
    [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(seqWght_prc(:,k,:),1);
    indSort_2 = hybrid(:,3);
    seqWght_prc_sort_2 = seqWght_prc(indSort_2,:,:);
    w_thres_sort = w_thres2(indSort_2);

    % take specific sequence weight
    selectSeqW = squeeze(seqWght_prc_sort_2(:,k,:));
    % take max weight and bin pos for each cell
    [selectSeqW_sum, selectSeqW_pos] = max(selectSeqW,[],2);
    % only consider putative principal cells
    prcInd = find(prcLabel_imec1 == 1); % get index value in original cell order
    prcID = clusterID_imec1(prcInd);
    depth_prc = depth_imec1(prcInd); % principle cell depth in original order
    prcIndSort = prcInd(indSort_2); % get the index value in seqNMF this sequence sorted order
    prcID_sort = prcID(indSort_2); % sorted cluster ID (ID is assigned in kilosort and phy)
    depth_prc_sort = depth_prc(indSort_2);  % principle cell depth in this sequence sorted order

    % we can do it in non-sorted original cell order way. But I want to see
    % sorted visualization
    prc_sort_Select_temp = find(selectSeqW_sum > w_thres_sort); % cell ind in sorted way
    % remove cells's position were the same at the sorting
    peakPos = selectSeqW_pos(prc_sort_Select_temp);
    remainInd = removeConsBin(peakPos,p.consBinThres); % remove bins if >= p.consBinThres in same bin

    prcID_sort_select_temp = prcID_sort(prc_sort_Select_temp); % cluster ID (ID is assigned in kilosort and phy)
    prcInd_sort_select_temp = prcIndSort(prc_sort_Select_temp); % cell Ind in all cells
    depth_Pyr_sort_select_temp = depth_prc_sort(prc_sort_Select_temp);
    pyrSelect_temp = indSort_2(prc_sort_Select_temp);
    selectSeqW_final_temp = selectSeqW_sum(prc_sort_Select_temp);

    seq_cellNumTemp = length(prc_sort_Select_temp);
    pyrInd_confirm_temp = zeros(seq_cellNumTemp,1);
    for m = 1:seq_cellNumTemp
        pyrInd_confirm_temp(m) = find(clusterID_imec1 == prcID_sort_select_temp(m)); % confirmation
    end
    if any(pyrInd_confirm_temp ~= prcInd_sort_select_temp)
        error('Sorting is wrong')
    end

    prc_sort_Select = prc_sort_Select_temp(remainInd);
    prcSelect = pyrSelect_temp(remainInd);
    prcInd_confirm = pyrInd_confirm_temp(remainInd);

    % confirmation
    seq_cellNum = length(prcInd_confirm);

    h = figure(1);
    h.Position = [100,100,900,900];
    subplot(2,1,1)
    WPlot(seqWght_pyr_sort(:,:,:),k);
    TITLE1 = sprintf('%d%s%d',p.L_time,'s-Seq-Delay-Pyr: ',k);
    TITLE2 = 'Sort by all seqs';
    title({TITLE1;TITLE2},'Interpreter','None')
    subplot(2,4,5)
    WPlot(seqWght_prc_sort_2(:,k,:),1);
    TITLE1 = sprintf('%d%s%d',p.L_time,'s-Seq: ',k);
    TITLE2 = 'Sort by all pyr cells';
    title({TITLE1;TITLE2},'Interpreter','None')
    subplot(2,4,6)
    plot(selectSeqW_sum,length(selectSeqW_sum):-1:1)
    hold on;
    plot(w_thres_sort,length(selectSeqW_sum):-1:1)

    subplot(2,4,7)
    WPlot(seqWght_prc_sort_2(prc_sort_Select_temp,k,:),1);
    hold on
    plot(peakPos(~remainInd),length(remainInd) - find(remainInd==0),'rx')
    TITLE1 = sprintf('%s%d%s%s%s','Seq: ',k);
    TITLE2 = sprintf('Activated cells, n=: %d',length(remainInd));
    title({TITLE1;TITLE2},'Interpreter','None')

    subplot(2,4,8)
    WPlot(seqWght_prc_sort_2(prc_sort_Select,k,:),1);
    TITLE1 = sprintf('%s%d%s%s%s','Seq: ',k);
    TITLE2 = sprintf('Final activated cells, n=: %d',length(prc_sort_Select));
    title({TITLE1;TITLE2},'Interpreter','None')


    % plot delay and spatial linear maps by selected cells in selected sequences
    % indSort_pyr
    h = figure(2);
    h.Position = [100,100,1600,900];
    subplot(2,5,1)
    WPlot(seqWght_prc_sort_2(prc_sort_Select,k,:),1);
    yl = ylim;
    TITLE1= sprintf('%d%s%d%s',p.L_time,'s-Seq-',k,'Selected prc Cells');
    title({TITLE1},'Interpreter','None')

    for j = 1:length(blockName2)
        subplot(2,5,j+1)
        timeMap_Pyr = timeMap_Def1.(blockName2{j});
        pre_Map_Pyr = pre_Map_Def1.(blockName2{j});
        post_Map_Pyr = post_Map_Def1.(blockName2{j});
        cellMapTemp = timeMap_Pyr(prcInd_confirm,:);
        pre_map1_sort = pre_Map_Pyr(prcInd_confirm,:);
        post_map1_sort = post_Map_Pyr(prcInd_confirm,:);
        cellMapTemp_Sort_Norm = [pre_map1_sort,cellMapTemp,post_map1_sort]./max(cellMapTemp,[],2);
        imagesc(flipud(cellMapTemp_Sort_Norm))
        colormap(jet)
        title(blockName2{j},'Interpreter','None')
        axis on
        set(gca, 'ydir', 'normal')
        set(gca, 'xtick', [1,size(pre_map1_sort,2),size([pre_map1_sort,cellMapTemp],2),size(cellMapTemp_Sort_Norm,2)]);
        if contains(blockName2{j},'10')
            set(gca, 'xticklabels', [-3 0 10 13]);
        else
            set(gca, 'xticklabels', [-3 0 30 33]);
        end
        clim([0 1])
        ylim(yl);
        xlabel('Time (s)')

    end

    subplot(2,5,6)
    linMapTemp = map_All_imec1;
    linMapTemp_Sort = linMapTemp(prcInd_confirm,:);
    linMapTemp_Sort_Norm_on = linMapTemp_Sort./max(linMapTemp_Sort,[],2);
    imagesc(flipud(linMapTemp_Sort_Norm_on))
    set(gca, 'ydir', 'normal')
    colormap(jet)
    set(gca, 'xtick', p.xTickLin);
    set(gca, 'xticklabels', p.xTickLabel);
    hold on; vertplot(p.linROI, 0, size(linMapTemp_Sort,1), 'r--');
    vertplot(p.linChoice, 0, size(linMapTemp_Sort,1), 'r--');
    TITLE1 = sprintf('%s','LinMap_All');
    title({TITLE1},'Interpreter','None')
    ylim(yl);

    if seq_cellNum >= p.min_pyrNum
        if any(prcInd_confirm ~= seqNMF_DelayDetect_Event_Extract.(seqName).cellInd)
            error('Sorting is wrong')
        end

        % plot ratster plot for each trial from selected seq NMF cells
        % look for best fit val and position

        % in the figure, on10_1 trial 2 and on30_2 trial 8 plotted
        blockInd = [1,7];
        trialInd = [8,2];
        orderVec = 1:seq_cellNum;
        h_thres_1 = seqNMF_DelayDetect_Event_Extract.(seqName).h_thres_1;

        for j = 1:length(blockInd)
            % plot all cells in the sequences in each trial
            % plot selected cells in the sequence in each trial
            % load analyzed positions
            load(fullfile(dataFolder,'Position',blockName{blockInd(j)},'ratesByECLR_imec1.mat'));
            delayZoneFile = fullfile(dataFolder,'Position',blockName{blockInd(j)}, 'Fig8DelayZonePos.mat');
            load(delayZoneFile);
            % stemZoneFile = fullfile(sessInfo(i).NIDQ,'Processed/Position',sessDirs{j}, 'Fig8StemZonePos.mat');
            % load(stemZoneFile);
            delayFile = fullfile(dataFolder,'Position',blockName{blockInd(j)}, 'PathZone.mat');
            load(delayFile);
            pathFile = fullfile(dataFolder,'Position',blockName{blockInd(j)}, 'pathData.mat');
            pathData = load(pathFile);
            pathData.t = posT_Sync_NIDQ.(blockName{blockInd(j)}); % synchronized neuralynx to neuropixels time
            pathLinFile = fullfile(dataFolder,'Position',blockName{blockInd(j)}, 'pathDataLinear.mat');
            pathLinData = load(pathLinFile);
            pathLinData.t = posT_Sync_NIDQ.(blockName{blockInd(j)});


            delay_Tstart = Fig8DelayZonePos.delayPos1.startT_Sync;
            stemStart = Fig8DelayZonePos.delayPos1.endT_Sync; % stem start ~= delay pos end
            returnStart = PathZone.posStartT_Sync.Return;
            baseStart = PathZone.posStartT_Sync.Base;
            baseEnd = PathZone.posEndT_Sync.Base;
            % choiceStart = Fig8StemZonePos.stemPos.endT; % choice start ~= stem pos end
            choiceEnd = PathZone.posEndT_Sync.Choice;
            rewardEnd = PathZone.posEndT_Sync.Reward;
            trialNum = length(delay_Tstart);


            if contains(blockName{blockInd(j)},'on')
                delaySpeed = Delay_Speed.(blockName{blockInd(j)}).speed;
            else
                delaySpeed = ones(1,trialNum);
            end

            if contains(blockName{blockInd(j)},'10')
                maxT = 10;
            else
                maxT = 30;
            end

            for m = trialInd(j) % for m = 1:trialNum    for all trials
                tridx = pathData.t >= returnStart(m) & pathData.t <= rewardEnd(m);
                posx = pathLinData.x(tridx);
                post = pathLinData.t(tridx);
                % find time and distance bins for return arm
                [~,posInd] = min(abs(pathLinData.t - baseEnd(m)));
                baseEndPos = pathLinData.x(posInd);
                bin_start.distance.return = [];
                bin_Width.distance.return = [];
                bin_start.time.return = [];
                bin_Width.time.return = [];
                for mm = 1:length(p.distanceBin_return)
                    binDistWidth = p.distanceBin_return(mm);
                    binDistTemp = 0:p.distIncre:(baseEndPos - binDistWidth);

                    for nn = 1:length(binDistTemp)
                        [diff_dist,ind_start] = min(abs(posx - binDistTemp(nn)));
                        [~,ind_end] = min(abs(posx - (binDistTemp(nn) + binDistWidth)));
                        if diff_dist <= p.distIncre
                            bin_start.time.return = [bin_start.time.return,post(ind_start)];
                            bin_Width.time.return = [bin_Width.time.return,post(ind_end) - post(ind_start)];

                            bin_start.distance.return = [bin_start.distance.return,binDistTemp(nn)];
                            bin_Width.distance.return = [bin_Width.distance.return,binDistWidth];
                        end
                    end
                end

                % find time and distance bins for stem arm
                stemStartT = stemStart(m);
                stemEndT = choiceEnd(m);
                [~,posInd] = min(abs(pathLinData.t - stemStartT));
                stemStartPos = pathLinData.x(posInd);
                [~,posInd] = min(abs(pathLinData.t - stemEndT));
                stemEndPos = pathLinData.x(posInd);

                bin_start.distance.stem = [];
                bin_Width.distance.stem = [];
                bin_start.time.stem = [];
                bin_Width.time.stem = [];

                for mm = 1:length(p.distanceBin_stem)
                    binDistWidth = p.distanceBin_stem(mm);
                    binDistTemp = stemStartPos:p.distIncre:(stemEndPos - binDistWidth);

                    for nn = 1:length(binDistTemp)
                        [diff_dist,ind_start] = min(abs(posx - binDistTemp(nn)));
                        [~,ind_end] = min(abs(posx - (binDistTemp(nn) + binDistWidth)));
                        if diff_dist <= p.distIncre
                            bin_start.time.stem = [bin_start.time.stem,post(ind_start)];
                            bin_Width.time.stem = [bin_Width.time.stem,post(ind_end) - post(ind_start)];

                            bin_start.distance.stem = [bin_start.distance.stem,binDistTemp(nn)];
                            bin_Width.distance.stem = [bin_Width.distance.stem,binDistWidth];
                        end
                    end
                end

                % left return arm
                h = figure;
                h.Position = [100,100,1600,800];

                % Top Panel
                % fig8 task
                subplot(2,3,1)
                tridx = pathData.t >= returnStart(m) & pathData.t <= rewardEnd(m);
                plot(pathData.x, pathData.y, 'color', ones(1,3)*.5);
                hold on
                plot(pathData.x(tridx), pathData.y(tridx), 'b');
                set (gca,'YDir','reverse')
                xlim([-100 100]); ylim([-100 100]); axis ij; axis square
                axis off
                TITLE1 = sprintf('%s-%s-Trial: %d%',seqName,blockName{j},m);
                title({TITLE1},'Interpreter','None');

                % this sequence
                subplot(2,3,3)
                WPlot(seqWght_prc_sort_2(prc_sort_Select,k,:),1);
                TITLE1 = sprintf('%d%s%s',p.L_time,'s-Seq: ',seqName);
                TITLE2 = sprintf('Final activated cells, n=: %d',length(prc_sort_Select));
                title({TITLE1;TITLE2},'Interpreter','None')

                h1 = subplot(2,3,2);
                seqH = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{blockInd(j)}).seqStrength_Delay{m};
                seqH_filter = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{blockInd(j)}).seqStrength_Delay_Newfit_filter{m};
                binTimeTemp = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{blockInd(j)}).filter_binTime{m};

                eventStartTs = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{blockInd(j)}).eventStartTs{m};
                eventEndTs = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{blockInd(j)}).eventEndTs{m};
                eventHPeak = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{blockInd(j)}).eventHPeak{m};

                plot((binTimeTemp+p.L_time/2-delay_Tstart(m)) * delaySpeed(m),seqH,'k');
                hold on
                plot((binTimeTemp+p.L_time/2-delay_Tstart(m)) * delaySpeed(m),seqH_filter,'b');
                plot((binTimeTemp([1,end])+p.L_time/2-delay_Tstart(m)) * delaySpeed(m),h_thres_1*[1,1],'--','Color',[1,0,0]);
                text(600,h_thres_1,sprintf('%1.4f',h_thres_1))

                if ~isempty(eventStartTs)
                    for kk = 1:length(eventStartTs)
                        [~,startInd] = min(abs(eventStartTs(kk) - binTimeTemp - p.L_time/2));
                        [~,endInd] = min(abs(eventEndTs(kk) - binTimeTemp - p.L_time/2));

                        plot((binTimeTemp(startInd:endInd)+ p.L_time/2-delay_Tstart(m)) * delaySpeed(m),seqH_filter(startInd:endInd),'r')
                    end
                end

                % return & stem
                for nn = 1:length(prcInd_confirm)
                    cellIndTemp = prcInd_confirm(nn);
                    tsp = pixelCluster_imec1.clusterInfo{cellIndTemp}.spkTs;
                    % return + base
                    tspInd = tsp >= returnStart(m) & tsp <= delay_Tstart(m);
                    tsp_plot = tsp(tspInd);
                    [spkx_return,spky,newTsp,spkPosInd] = spk2posInd(tsp_plot,pathLinData.x,pathLinData.y,pathLinData.t);

                    subplot(2,3,4)
                    % take only ascending distance spikes to avoid animal run back spikes
                    spkx_return_ascend = spkx_return;
                    removeInd = zeros(size(spkx_return));

                    for indtemp = 2:length(spkx_return)
                        if spkx_return_ascend(indtemp) < spkx_return_ascend(indtemp-1)
                            spkx_return_ascend(indtemp) = spkx_return_ascend(indtemp-1);
                            removeInd(indtemp) = 1;
                        end
                    end
                    spkx_return_ascend(removeInd==1) = [];

                    xPoints = [spkx_return_ascend';spkx_return_ascend'];
                    ypos = nn;
                    yPoints = [ypos+zeros(size(spkx_return_ascend'))-0.3;ypos+zeros(size(spkx_return_ascend'))+0.3];


                    if ~isempty(tsp_plot)
                        plot(xPoints,yPoints,'k')
                    end
                    hold on


                    % delay
                    tspInd = tsp >= delay_Tstart(m) & tsp <= stemStart(m);
                    tsp_plot = tsp(tspInd) - delay_Tstart(m) ;
                    spkx_delay = tsp_plot * delaySpeed(m);
                    h2 = subplot(2,3,5);
                    xPoints = [spkx_delay';spkx_delay'];
                    ypos = nn;
                    yPoints = [ypos+zeros(size(tsp_plot'))-0.3;ypos+zeros(size(tsp_plot'))+0.3];
                    if ~isempty(tsp_plot)
                        plot(xPoints,yPoints,'k')
                    end
                    hold on

                    % stem + choice + reward
                    tspInd = tsp >= stemStart(m) & tsp <= rewardEnd(m);
                    tsp_plot = tsp(tspInd);
                    [spkx_stem,spky,newTsp,spkPosInd] = spk2posInd(tsp_plot,pathLinData.x,pathLinData.y,pathLinData.t);
                    subplot(2,3,6)
                    % take only ascending distance spikes to avoid animal run back spikes
                    spkx_stem_ascend = spkx_stem;
                    removeInd = zeros(size(spkx_stem));
                    for indtemp = 2:length(spkx_stem)
                        if spkx_stem_ascend(indtemp) < spkx_stem_ascend(indtemp-1)
                            spkx_stem_ascend(indtemp) = spkx_stem_ascend(indtemp-1);
                            removeInd(indtemp) = 1;
                        end
                    end
                    spkx_stem_ascend(removeInd==1) = [];

                    xPoints = [spkx_stem_ascend';spkx_stem_ascend'];
                    ypos = nn;
                    yPoints = [ypos+zeros(size(spkx_stem_ascend'))-0.3;ypos+zeros(size(spkx_stem_ascend'))+0.3];

                    if ~isempty(tsp_plot)
                        plot(xPoints,yPoints,'k')
                    end
                    hold on
                end

                % delay events
                subplot(2,3,5)
                if ~isempty(eventStartTs)
                    for kk = 1:length(eventStartTs)
                        for nn = 1:length(prcInd_confirm)
                            cellInd = prcInd_confirm(nn);
                            tsp = pixelCluster_imec1.clusterInfo{cellInd}.spkTs;
                            tspInd = tsp >= eventStartTs(kk) & tsp <= eventEndTs(kk);

                            tsp_plot = tsp(tspInd) - delay_Tstart(m) ;
                            spkx_delay = tsp_plot * delaySpeed(m);
                            xPoints = [spkx_delay';spkx_delay'];
                            ypos = nn;
                            yPoints = [ypos+zeros(size(tsp_plot'))-0.3;ypos+zeros(size(tsp_plot'))+0.3];
                            if ~isempty(tsp_plot)
                                plot(xPoints,yPoints,'Color',[1,0.4,0])
                            end
                            hold on
                        end
                    end
                end
                xline(maxT*delaySpeed(m),'r--');
                linkaxes([h1,h2],'x');

                % return and stem
                % calculate the fit, go through every dist bin and calculate
                % 1. spike count 2. spike median time 3. spike mean time 4. overall order
                %                 bin_start.distance.return = [];
                %                 bin_Width.distance.return = [];
                %                 bin_start.time.return = [];
                %                 bin_Width.time.return = [];
                spike_Count_return = zeros(length(bin_start.time.return),length(prcInd_confirm));
                spike_ts_Median_return = zeros(length(bin_start.time.return),length(prcInd_confirm));
                spike_ts_Mean_return = zeros(length(bin_start.time.return),length(prcInd_confirm));
                spike_ts_MedianOrder_return = zeros(length(bin_start.time.return),length(prcInd_confirm));
                spike_ts_MeanOrder_return = zeros(length(bin_start.time.return),length(prcInd_confirm));
                spike_ts_Corr_rank_return = nan(length(bin_start.time.return),1);
                spike_ts_Corr_rank_p_return = nan(length(bin_start.time.return),1);

                spike_dist_Median_return = zeros(length(bin_start.time.return),length(prcInd_confirm));
                spike_dist_Mean_return = zeros(length(bin_start.time.return),length(prcInd_confirm));
                spike_dist_MedianOrder_return = zeros(length(bin_start.time.return),length(prcInd_confirm));
                spike_dist_MeanOrder_return = zeros(length(bin_start.time.return),length(prcInd_confirm));
                spike_dist_Corr_rank_return = nan(length(bin_start.time.return),1);
                spike_dist_Corr_rank_p_return = nan(length(bin_start.time.return),1);
                spike_dist_Corr_rank_weighted_return = nan(length(bin_start.time.return),1);

                for mm = 1:length(bin_start.time.return)
                    startT = bin_start.time.return(mm);
                    endT = bin_start.time.return(mm) + bin_Width.time.return(mm);

                    spike_ts_CountTemp = zeros(1,length(prcInd_confirm));
                    spike_ts_MedianTemp = nan(1,length(prcInd_confirm));
                    spike_ts_MeanTemp = nan(1,length(prcInd_confirm));

                    spike_dist_MedianTemp = nan(1,length(prcInd_confirm));
                    spike_dist_MeanTemp = nan(1,length(prcInd_confirm));

                    for nn = 1:length(prcInd_confirm)
                        cellInd = prcInd_confirm(nn);
                        tsp = pixelCluster_imec1.clusterInfo{cellInd}.spkTs;
                        tspInd = tsp >= startT & tsp <= endT;
                        % spike time
                        spike_ts_CountTemp(nn) = sum(tspInd);
                        spike_ts_MedianTemp(nn) = median(tsp(tspInd));
                        spike_ts_MeanTemp(nn) = nanmean(tsp(tspInd));

                        % spike distance
                        [spkx_return,~,~,~] = spk2posInd(tsp(tspInd),pathLinData.x,pathLinData.y,pathLinData.t);
                        spkx_return_ascend = spkx_return;
                        removeInd = zeros(size(spkx_return));

                        for indtemp = 2:length(spkx_return)
                            if spkx_return_ascend(indtemp) < spkx_return_ascend(indtemp-1)
                                spkx_return_ascend(indtemp) = spkx_return_ascend(indtemp-1);
                                removeInd(indtemp) = 1;
                            end
                        end
                        spkx_return_ascend(removeInd==1) = [];

                        spike_dist_MedianTemp(nn) = nanmedian(spkx_return_ascend);
                        spike_dist_MeanTemp(nn) = nanmean(spkx_return_ascend);
                    end

                    nanInd = spike_ts_CountTemp == 0;
                    % time
                    [~,order_medianTemp] = sort(spike_ts_MedianTemp(~nanInd));
                    rankTemp = zeros(1,sum(~nanInd));
                    rankTemp(order_medianTemp) = 1:numel(order_medianTemp);
                    order_ts_median = nan(size(spike_ts_MedianTemp));
                    order_ts_median(~nanInd) = rankTemp;

                    [~,order_meanTemp] = sort(spike_ts_MeanTemp(~nanInd));
                    rankTemp = zeros(1,sum(~nanInd));
                    rankTemp(order_meanTemp) = 1:numel(order_meanTemp);
                    order_ts_mean = nan(size(spike_ts_MeanTemp));
                    order_ts_mean(~nanInd) = rankTemp;

                    spike_Count_return(mm,:) = spike_ts_CountTemp;
                    spike_ts_Median_return(mm,:) = spike_ts_MedianTemp;
                    spike_ts_Mean_return(mm,:) = spike_ts_MeanTemp;

                    spike_ts_MedianOrder_return(mm,:) = order_ts_median;
                    spike_ts_MeanOrder_return(mm,:) = order_ts_mean;

                    % distance
                    [~,order_medianTemp] = sort(spike_dist_MedianTemp(~nanInd));
                    rankTemp = zeros(1,sum(~nanInd));
                    rankTemp(order_medianTemp) = 1:numel(order_medianTemp);
                    order_dist_median = nan(size(spike_dist_MedianTemp));
                    order_dist_median(~nanInd) = rankTemp;

                    [~,order_meanTemp] = sort(spike_dist_MeanTemp(~nanInd));
                    rankTemp = zeros(1,sum(~nanInd));
                    rankTemp(order_meanTemp) = 1:numel(order_meanTemp);
                    order_dist_mean = nan(size(spike_dist_MeanTemp));
                    order_dist_mean(~nanInd) = rankTemp;

                    spike_dist_Median_return(mm,:) = spike_dist_MedianTemp;
                    spike_dist_Mean_return(mm,:) = spike_dist_MeanTemp;

                    spike_dist_MedianOrder_return(mm,:) = order_dist_median;
                    spike_dist_MeanOrder_return(mm,:) = order_dist_mean;

                    % calculated the correlation rank order
                    if sum(~isnan(order_ts_median))>= p.minCorrCell

                        % time
                        [spike_ts_Corr_rank_return(mm),spike_ts_Corr_rank_p_return(mm)] = corr(orderVec',order_ts_mean', 'Type', 'Spearman', 'Rows', 'complete');

                        % distance
                        [spike_dist_Corr_rank_return(mm),spike_dist_Corr_rank_p_return(mm)] = corr(orderVec',order_dist_mean', 'Type', 'Spearman', 'Rows', 'complete');

                    end
                end

                spikeCorr_rank_validIdx = find(spike_dist_Corr_rank_p_return <= 0.05);
                spikeCorr_rank_Valid = spike_dist_Corr_rank_return(spikeCorr_rank_validIdx);
                spikeCorr_rank_p_Valid = spike_dist_Corr_rank_p_return(spikeCorr_rank_validIdx);

                [maxCorr,maxInd] = max(spikeCorr_rank_Valid);
                maxCorr_p = spikeCorr_rank_p_Valid(maxInd);

                if maxCorr > 0
                    maxStartT = bin_start.time.return(spikeCorr_rank_validIdx(maxInd));
                    maxEndT = bin_start.time.return(spikeCorr_rank_validIdx(maxInd))+bin_Width.time.return(spikeCorr_rank_validIdx(maxInd));
                else
                    maxStartT = [];
                    maxEndT = [];
                end

                [minCorr,minInd] = min(spikeCorr_rank_Valid);
                minCorr_p =spikeCorr_rank_p_Valid(minInd);
                if minCorr < 0
                    minStartT = bin_start.time.return(spikeCorr_rank_validIdx(minInd));
                    minEndT = bin_start.time.return(spikeCorr_rank_validIdx(minInd))+bin_Width.time.return(spikeCorr_rank_validIdx(minInd));
                else
                    minStartT = [];
                    minEndT = [];
                end


                subplot(2,3,4)
                for nn = 1:length(prcInd_confirm)
                    cellIndTemp = prcInd_confirm(nn);
                    tsp = pixelCluster_imec1.clusterInfo{cellIndTemp}.spkTs;
                    % % max corr
                    % if ~isempty(maxStartT)
                    %     tspInd = tsp >= maxStartT & tsp <= maxEndT;
                    %     tsp_plot = tsp(tspInd);
                    %     [spkx_return,spky,newTsp,spkPosInd] = spk2posInd(tsp_plot,pathLinData.x,pathLinData.y,pathLinData.t);
                    %     xPoints = [spkx_return';spkx_return'];
                    %     ypos = nn;
                    %     yPoints = [ypos+zeros(size(spkx_return'))-0.3;ypos+zeros(size(spkx_return'))+0.3];
                    %     if ~isempty(tsp_plot)
                    %         plot(xPoints,yPoints,'Color',[1,0.4,0])
                    %     end
                    % end

                    % min corr
                    if ~isempty(minStartT)
                        % return + base
                        tspInd = tsp >= returnStart(m) & tsp <= delay_Tstart(m);
                        tsp_plot = tsp(tspInd);
                        [spkx_return,spky,newTsp,spkPosInd] = spk2posInd(tsp_plot,pathLinData.x,pathLinData.y,pathLinData.t);

                        subplot(2,3,4)
                        % take only ascending distance spikes to avoid animal run back spikes
                        spkx_return_ascend = spkx_return;
                        removeInd = zeros(size(spkx_return));

                        for indtemp = 2:length(spkx_return)
                            if spkx_return_ascend(indtemp) < spkx_return_ascend(indtemp-1)
                                spkx_return_ascend(indtemp) = spkx_return_ascend(indtemp-1);
                                removeInd(indtemp) = 1;
                            end
                        end
                        spkx_return_ascend(removeInd==1) = [];
                        newTsp(removeInd==1) = [];

                        newInd = newTsp >= minStartT & newTsp <= minEndT;
                        spkx_return_plot = spkx_return_ascend(newInd);

                        xPoints = [spkx_return_plot';spkx_return_plot'];
                        ypos = nn;
                        yPoints = [ypos+zeros(size(spkx_return_plot'))-0.3;ypos+zeros(size(spkx_return_plot'))+0.3];
                        if ~isempty(tsp_plot)
                            plot(xPoints,yPoints,'Color',[1,0.4,0])
                        end
                    end
                end

                % if ~isempty(maxStartT)
                %     TEXT = sprintf('max Corr: %1.2f',maxCorr);
                %     text(bin_start.distance.return(maxInd),0,TEXT)
                %     TEXT = sprintf('pVal: %.3e',maxCorr_p);
                %     text(bin_start.distance.return(maxInd),2,TEXT)
                % end

                if ~isempty(minStartT)
                    TEXT = sprintf('min Corr: %1.2f',minCorr);
                    text(bin_start.distance.return(minInd),0,TEXT)
                    TEXT = sprintf('pVal: %.3e',minCorr_p);
                    text(bin_start.distance.return(minInd),2,TEXT)
                end

                % stem
                spike_Count_stem = zeros(length(bin_start.time.stem),length(prcInd_confirm));
                spike_ts_Median_stem = zeros(length(bin_start.time.stem),length(prcInd_confirm));
                spike_ts_Mean_stem = zeros(length(bin_start.time.stem),length(prcInd_confirm));
                spike_ts_MedianOrder_stem = zeros(length(bin_start.time.stem),length(prcInd_confirm));
                spike_ts_MeanOrder_stem = zeros(length(bin_start.time.stem),length(prcInd_confirm));
                spike_ts_Corr_rank_stem = nan(length(bin_start.time.stem),1);
                spike_ts_Corr_rank_p_stem = nan(length(bin_start.time.stem),1);

                spike_dist_Median_stem = zeros(length(bin_start.time.stem),length(prcInd_confirm));
                spike_dist_Mean_stem = zeros(length(bin_start.time.stem),length(prcInd_confirm));
                spike_dist_MedianOrder_stem = zeros(length(bin_start.time.stem),length(prcInd_confirm));
                spike_dist_MeanOrder_stem = zeros(length(bin_start.time.stem),length(prcInd_confirm));
                spike_dist_Corr_rank_stem = nan(length(bin_start.time.stem),1);
                spike_dist_Corr_rank_p_stem = nan(length(bin_start.time.stem),1);
                spike_dist_Corr_rank_weighted_stem = nan(length(bin_start.time.stem),1);

                for mm = 1:length(bin_start.time.stem)
                    startT = bin_start.time.stem(mm);
                    endT = bin_start.time.stem(mm) + bin_Width.time.stem(mm);

                    spike_ts_CountTemp = zeros(1,length(prcInd_confirm));
                    spike_ts_MedianTemp = nan(1,length(prcInd_confirm));
                    spike_ts_MeanTemp = nan(1,length(prcInd_confirm));

                    spike_dist_MedianTemp = nan(1,length(prcInd_confirm));
                    spike_dist_MeanTemp = nan(1,length(prcInd_confirm));

                    for nn = 1:length(prcInd_confirm)
                        cellInd = prcInd_confirm(nn);
                        tsp = pixelCluster_imec1.clusterInfo{cellInd}.spkTs;
                        tspInd = tsp >= startT & tsp <= endT;
                        % spike time
                        spike_ts_CountTemp(nn) = sum(tspInd);
                        spike_ts_MedianTemp(nn) = median(tsp(tspInd));
                        spike_ts_MeanTemp(nn) = nanmean(tsp(tspInd));

                        % spike distance
                        [spkx_stem,~,~,~] = spk2posInd(tsp(tspInd),pathLinData.x,pathLinData.y,pathLinData.t);
                        % take only ascending distance spikes to avoid animal run back spikes
                        spkx_stem_ascend = spkx_stem;
                        removeInd = zeros(size(spkx_stem));
                        for indtemp = 2:length(spkx_stem)
                            if spkx_stem_ascend(indtemp) < spkx_stem_ascend(indtemp-1)
                                spkx_stem_ascend(indtemp) = spkx_stem_ascend(indtemp-1);
                                removeInd(indtemp) = 1;
                            end
                        end
                        spkx_stem_ascend(removeInd==1) = [];

                        spike_dist_MedianTemp(nn) = nanmedian(spkx_stem_ascend);
                        spike_dist_MeanTemp(nn) = nanmean(spkx_stem_ascend);
                    end

                    nanInd = spike_ts_CountTemp == 0;
                    % time
                    [~,order_medianTemp] = sort(spike_ts_MedianTemp(~nanInd));
                    rankTemp = zeros(1,sum(~nanInd));
                    rankTemp(order_medianTemp) = 1:numel(order_medianTemp);
                    order_ts_median = nan(size(spike_ts_MedianTemp));
                    order_ts_median(~nanInd) = rankTemp;

                    [~,order_meanTemp] = sort(spike_ts_MeanTemp(~nanInd));
                    rankTemp = zeros(1,sum(~nanInd));
                    rankTemp(order_meanTemp) = 1:numel(order_meanTemp);
                    order_ts_mean = nan(size(spike_ts_MeanTemp));
                    order_ts_mean(~nanInd) = rankTemp;

                    spike_Count_stem(mm,:) = spike_ts_CountTemp;
                    spike_ts_Median_stem(mm,:) = spike_ts_MedianTemp;
                    spike_ts_Mean_stem(mm,:) = spike_ts_MeanTemp;

                    spike_ts_MedianOrder_stem(mm,:) = order_ts_median;
                    spike_ts_MeanOrder_stem(mm,:) = order_ts_mean;

                    % distance
                    [~,order_medianTemp] = sort(spike_dist_MedianTemp(~nanInd));
                    rankTemp = zeros(1,sum(~nanInd));
                    rankTemp(order_medianTemp) = 1:numel(order_medianTemp);
                    order_dist_median = nan(size(spike_dist_MedianTemp));
                    order_dist_median(~nanInd) = rankTemp;

                    [~,order_meanTemp] = sort(spike_dist_MeanTemp(~nanInd));
                    rankTemp = zeros(1,sum(~nanInd));
                    rankTemp(order_meanTemp) = 1:numel(order_meanTemp);
                    order_dist_mean = nan(size(spike_dist_MeanTemp));
                    order_dist_mean(~nanInd) = rankTemp;

                    spike_dist_Median_stem(mm,:) = spike_dist_MedianTemp;
                    spike_dist_Mean_stem(mm,:) = spike_dist_MeanTemp;

                    spike_dist_MedianOrder_stem(mm,:) = order_dist_median;
                    spike_dist_MeanOrder_stem(mm,:) = order_dist_mean;

                    % calculated the correlation
                    % rank order 
                    if sum(~isnan(order_ts_median))>= p.minCorrCell

                        % time
                        [spike_ts_Corr_rank_stem(mm),spike_ts_Corr_rank_p_stem(mm)] = corr(orderVec',order_ts_mean', 'Type', 'Spearman', 'Rows', 'complete');

                        % distance
                        [spike_dist_Corr_rank_stem(mm),spike_dist_Corr_rank_p_stem(mm)] = corr(orderVec',order_dist_mean', 'Type', 'Spearman', 'Rows', 'complete');

                    end
                end

                spikeCorr_rank_validIdx = find(spike_dist_Corr_rank_p_stem <= 0.05);
                spikeCorr_rank_Valid = spike_dist_Corr_rank_stem(spikeCorr_rank_validIdx);
                spikeCorr_rank_p_Valid = spike_dist_Corr_rank_p_stem(spikeCorr_rank_validIdx);

                [maxCorr,maxInd] = max(spikeCorr_rank_Valid);
                maxCorr_p = spikeCorr_rank_p_Valid(maxInd);

                if maxCorr > 0
                    maxStartT = bin_start.time.stem(spikeCorr_rank_validIdx(maxInd));
                    maxEndT = bin_start.time.stem(spikeCorr_rank_validIdx(maxInd))+bin_Width.time.stem(spikeCorr_rank_validIdx(maxInd));
                else
                    maxStartT = [];
                    maxEndT = [];
                end

                [minCorr,minInd] = min(spikeCorr_rank_Valid);
                minCorr_p = spikeCorr_rank_p_Valid(minInd);
                if minCorr < 0
                    minStartT = bin_start.time.stem(spikeCorr_rank_validIdx(minInd));
                    minEndT = bin_start.time.stem(spikeCorr_rank_validIdx(minInd))+bin_Width.time.stem(spikeCorr_rank_validIdx(minInd));
                else
                    minStartT = [];
                    minEndT = [];
                end


                subplot(2,3,6)
                for nn = 1:length(prcInd_confirm)
                    cellIndTemp = prcInd_confirm(nn);
                    tsp = pixelCluster_imec1.clusterInfo{cellIndTemp}.spkTs;
                    % max corr
                    if ~isempty(maxStartT)

                        tspInd = tsp >= stemStart(m) & tsp <= rewardEnd(m);
                        tsp_plot = tsp(tspInd);
                        [spkx_stem,spky,newTsp,spkPosInd] = spk2posInd(tsp_plot,pathLinData.x,pathLinData.y,pathLinData.t);

                        % take only ascending distance spikes to avoid animal run back spikes
                        spkx_stem_ascend = spkx_stem;
                        removeInd = zeros(size(spkx_stem));
                        for indtemp = 2:length(spkx_stem)
                            if spkx_stem_ascend(indtemp) < spkx_stem_ascend(indtemp-1)
                                spkx_stem_ascend(indtemp) = spkx_stem_ascend(indtemp-1);
                                removeInd(indtemp) = 1;
                            end
                        end
                        spkx_stem_ascend(removeInd==1) = [];
                        newTsp(removeInd==1) = [];

                        subplot(2,3,6)
                        newInd = newTsp >= maxStartT & newTsp <= maxEndT;
                        spkx_stem_plot = spkx_stem_ascend(newInd);

                        xPoints = [spkx_stem_plot';spkx_stem_plot'];
                        ypos = nn;
                        yPoints = [ypos+zeros(size(spkx_stem_plot'))-0.3;ypos+zeros(size(spkx_stem_plot'))+0.3];
                        if ~isempty(tsp_plot)
                            plot(xPoints,yPoints,'Color',[1,0.4,0])
                        end
                    end

                    % % min corr
                    % if ~isempty(minStartT)
                    %     tspInd = tsp >= minStartT & tsp <= minEndT;
                    %     tsp_plot = tsp(tspInd);
                    %     [spkx_return,spky,newTsp,spkPosInd] = spk2posInd(tsp_plot,pathLinData.x,pathLinData.y,pathLinData.t);
                    %     xPoints = [spkx_return';spkx_return'];
                    %     ypos = nn;
                    %     yPoints = [ypos+zeros(size(spkx_return'))-0.3;ypos+zeros(size(spkx_return'))+0.3];
                    %     if ~isempty(tsp_plot)
                    %         plot(xPoints,yPoints,'Color',[1,0.4,0])
                    %     end
                    % end
                end

                if ~isempty(maxStartT)
                    TEXT = sprintf('max Corr: %1.2f',maxCorr);
                    text(bin_start.distance.stem(maxInd),0,TEXT)
                    TEXT = sprintf('pVal: %.3e',maxCorr_p);
                    text(bin_start.distance.stem(maxInd),2,TEXT)
                end

                % if ~isempty(minStartT)
                %     TEXT = sprintf('min Corr: %1.2f',minCorr);
                %     text(bin_start.distance.return(minInd),0,TEXT)
                %     TEXT = sprintf('pVal: %.3e',minCorr_p);
                %     text(bin_start.distance.return(minInd),2,TEXT)
                % end



                subplot(2,3,4)
                ylim([0 length(prcInd_confirm)+1])
                set(gca, 'YTick', [1,length(prcInd_confirm)], 'YTickLabel', num2cell([1,length(prcInd_confirm)]), 'TickLength', [0, 0]);
                set(gca,'YDir','reverse');
                XLABEL = sprintf('Return + base (cm), %3.1f s', delay_Tstart(m)-returnStart(m));
                xlabel(XLABEL);
                TITLE1 = sprintf('%s-Trial: %d return',blockName{j},m);
                title(TITLE1);
                axis tight
                xlim([0 200])

                subplot(2,3,5)
                % xlim([0 190])
                ylim([0 length(prcInd_confirm)+1])
                set(gca, 'YTick', [1,length(prcInd_confirm)], 'YTickLabel', num2cell([1,length(prcInd_confirm)]), 'TickLength', [0, 0]);
                set(gca,'YDir','reverse');
                XLABEL = sprintf('Delay total (cm), %3.1f s', stemStart(m)-delay_Tstart(m));
                xlabel(XLABEL);
                TITLE1 = 'Delay';
                title(TITLE1);
                axis tight
                xlim([0 600])

                subplot(2,3,6)
                ylim([0 length(prcInd_confirm)+1])
                set(gca, 'YTick', [1,length(prcInd_confirm)], 'YTickLabel', num2cell([1,length(prcInd_confirm)]), 'TickLength', [0, 0]);
                set(gca,'YDir','reverse');
                XLABEL = sprintf('Stem + choice + reward (cm), %3.1f s', rewardEnd(m)-stemStart(m));
                xlabel(XLABEL);
                axis tight
                xlim([230 430])

            end
        end
    end
end
fprintf('Finished figure 6\n');


