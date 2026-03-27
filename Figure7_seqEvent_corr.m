% this code calculate and plot the auto- and cross-correlogram of sequence
% events
% Li Yuan, UCSD

clear all
close all

codePath = 'C:\Users\Li\Documents\MATLAB\DualProbe\CellProperty\Local\MEC_manuscript_code'; % change the folder to your local path
dataFolder = fullfile(codePath,'Data','1138_20240404');
addpath(genpath(codePath));

p.timeBin = 200/10^3;
p.L_time = 3; % unit in second
p.seqThresPct = 99;
p.consBinThres = 6;
p.w_min = 1;
p.min_pyrNum = 10; % here make it low to include more sequences, will increase to 15 for quant

p.savePlot = 1;
p.writeToFile = 0;

% 8 blocks in a recording session
blockName = {'on10_1','off10_1','on30_1','off30_1','on10_2','off10_2','on30_2','off30_2'};

% load synchrinized pos t
load(fullfile(dataFolder, 'posT_Sync_NIDQ.mat'));
% load cell int pyr label in MEC
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


% get valid cell ind
%     pyrLabel_imec1 = mec_clusterType.labelInd;
prcLabel_imec1 = Fig8_seq_Delay_detect.imec1.pyrInd;
avgRate_imec1 = mec_clusterType.avgRate;
depth_imec1 = mec_clusterType.clusterDepth;
cellLayer_imec1 = mec_clusterType.layerInd;
depthLayer_imec1 = mec_clusterType.posInd;
clusterID_imec1 = mec_clusterType.clusterID;
clusterNum = mec_clusterType.clusterNum;

% get weight dist for all pyr cells
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

validSeqCount = 0;
validSeq = [];
validWeight = [];
validPrcInd = [];

for k = 1:seqNum
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
    seq_cellNum = length(prcInd_confirm);

    if seq_cellNum >= p.min_pyrNum
        if any(prcInd_confirm ~= seqNMF_DelayDetect_Event_Extract.(seqName).cellInd)
            error('Sorting is wrong')
        else
            validSeqCount = validSeqCount + 1;
            validSeq = [validSeq,k];
            validWeight{validSeqCount} = seqWght_prc_sort_2(prc_sort_Select,k,:);
            validPrcInd{validSeqCount} = prcInd_confirm;
        end
    end
end

% to avoid delay concatenations
binNumInsert = ceil(30./p.timeBin);
binInsert = zeros(1,binNumInsert);

%% calculate and plot the correlogram
if validSeqCount >= 1
    seqBinH = cell(validSeqCount,1);
    seqBinH_threshold = cell(validSeqCount,1);
    seqTime = cell(validSeqCount,1);
    seqEvent_H = cell(validSeqCount,1);
    % get each sequence event time
    for j = [1,3,5,7]

        delayZoneFile = fullfile(dataFolder,'Position',blockName{j}, 'Fig8DelayZonePos.mat');
        load(delayZoneFile);
        delay_Tstart = Fig8DelayZonePos.delayPos1.startT_Sync;
        trialNum = length(delay_Tstart);

        for m = 1:trialNum

            % this sequence

            for k = 1:validSeqCount
                validInd = validSeq(k);
                seqName = sprintf('seq%d',validInd);
                % seqWeight = validWeight{k};
                % pyrInd_confirm = seqNMF_DelayDetect_Event_Extract.(seqName).cellInd;

                delayH = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{j}).seqStrength_Delay_Newfit_filter{m};
                threshold = seqNMF_DelayDetect_Event_Extract.(seqName).h_thres_1;
                delayH_threshold = delayH/threshold;
                eventStartTs = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{j}).eventStartTs{m};
                eventEndTs = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{j}).eventEndTs{m};
                eventPeakTs = (eventStartTs + eventEndTs)/2;
                eventHPeak = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{j}).eventHPeak{m};

                seqBinH{k} = [seqBinH{k},binInsert,delayH];
                seqBinH_threshold{k} = [seqBinH_threshold{k},binInsert,delayH_threshold];
                seqTime{k} = [seqTime{k},eventPeakTs];
                seqEvent_H{k} = [seqEvent_H{k},eventHPeak];

            end

        end
    end

    h1 = figure(1);
    h1.Position = [100,100,900,900];

    for k = 1:4
        validInd_1 = validSeq(k);
        seqName_1 = sprintf('seq%d',validInd_1);
        seqWeight_1 = validWeight{k};
        prcInd_confirm = seqNMF_DelayDetect_Event_Extract.(seqName_1).cellInd;

        figure(1)
        subplot(5,5,k+1)
        WPlot(seqWeight_1,1);
        TITLE1 = sprintf('%d%s%s',p.L_time,'s-',seqName_1);
        title({TITLE1},'Interpreter','None')

        subplot(5,5,5*k+1)
        WPlot(seqWeight_1,1);
        TITLE1 = sprintf('%d%s%s',p.L_time,'s-',seqName_1);
        TITLE2 = sprintf('%d cells',length(prcInd_confirm));
        title({TITLE1;TITLE2},'Interpreter','None')

    end

    maxLag = ceil(30 / p.timeBin);

    pairDiff = cell(validSeqCount);
    pairDiff_Norm = cell(validSeqCount);
    edge_pairs = cell(validSeqCount);
    p_val_pairs = nan(validSeqCount);

    for k = 1:4
        validInd_1 = validSeq(k);
        seqName_1 = sprintf('seq%d',validInd_1);
        seqWeight_1 = validWeight{k};

        seqH_1 = seqBinH{k};
        seqH_norm_1 = seqBinH_threshold{k};
        seqTime_1 = seqTime{k};

        for m = 1:k
            validInd_2 = validSeq(k);
            seqName = sprintf('seq%d',validInd_2);
            seqWeight_2 = validWeight{k};

            seqH_2 = seqBinH{m};
            seqH_norm_2 = seqBinH_threshold{m};
            seqTime_2 = seqTime{m};

            figure(1)
            subplot(5,5,k*5+m+1)
            % Pairwise differences (t2 - t1)
            d = bsxfun(@minus, seqTime_2(:), seqTime_1(:)');  % matrix of differences
            d = d(:);                                         % vector of differences

            % Lag bins
            lagWidth = 1;                  % bin width in seconds
            lagRange = -30:lagWidth:30;     % edges (choose symmetric range)
            [counts, edges] = histcounts(d, [lagRange-lagWidth/2, lagRange(end)+lagWidth/2]);
            lags = edges(1:end-1) + lagWidth/2;
            edge_pairs{k,m} = lags;

            % Optionally normalize to rate or number of reference events
            counts_norm = counts / (numel(seqTime_1) * lagWidth);  % rate per second per reference event


            % test significance
            centerWin = [ceil(-30/lagWidth) ceil(30/lagWidth)];   % change if your timescale is different
            centerIdx = lags >= centerWin(1) & lags <= centerWin(2);

            % simple metric: peak height in that window
            peak_real = max(counts_norm(centerIdx));

            % optional: baseline-subtracted peak (using flanks)
            flankIdx = lags < centerWin(1) | lags > centerWin(2);
            baseline = mean(counts_norm(flankIdx));
            metric_real = peak_real;  % this is your "sequence relatedness" metric

            T = max([seqTime_1(:); seqTime_2(:)]);
            nShuf = 10000;

            metric_shuf = zeros(nShuf,1);

            for s = 1:nShuf
                shift = rand * 15;                                % random shift in [0, T)
                seqTime_2_sh = mod(seqTime_2 + shift, T);        % circular wrap

                % recompute correlogram for shuffled times
                d_sh = bsxfun(@minus, seqTime_2_sh(:), seqTime_1(:)');
                d_sh = d_sh(:);

                [counts_sh, ~] = histcounts(d_sh, [lagRange-lagWidth/2, lagRange(end)+lagWidth/2]);
                counts_sh_norm = counts_sh / (numel(seqTime_1) * lagWidth);

                peak_sh = max(counts_sh_norm(centerIdx));
                baseline_sh = mean(counts_sh_norm(flankIdx));
                metric_shuf(s) = peak_sh;
            end

            p_val = (1 + sum(metric_shuf > metric_real)) / (nShuf + 1);

            pairDiff{k,m} = counts;
            pairDiff_Norm{k,m} = counts_norm;
            p_val_pairs(k,m) = p_val;


            % Plot
            bar(lags, counts_norm)
            xlabel('Lag (s)')
            ylabel('Counts')
            title('pairWise event')
            TEXT = sprintf('%.2e',p_val);
            text(10,0.2,num2str(TEXT))
        end
    end
end

%% plot single trial delay example
% in the figure, on10_1 trial 2 and on30_2 trial 8 plotted
blockInd = [7,7];
trialInd = [5,8];
colorVal = [0,0.2,1;1,0.6,0.2;0.4,0.2,0.8;0,0.4,0];
for j = 1:length(blockInd)
    h = figure;
    h.Position = [100,100,1600,800];

    % plot all cells in the sequences in each trial
    % plot selected cells in the sequence in each trial
    % load analyzed positions
    delayZoneFile = fullfile(dataFolder,'Position',blockName{blockInd(j)}, 'Fig8DelayZonePos.mat');
    load(delayZoneFile);

    delay_Tstart = Fig8DelayZonePos.delayPos1.startT_Sync;
    stemStart = Fig8DelayZonePos.delayPos1.endT_Sync;
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

        for k = 1:4
            seqName = sprintf('seq%d',k);
            prcInd_confirm = validPrcInd{k};

            % this sequence
            subplot(4,4,4*(k-1)+1)
            seqWeight_1 = validWeight{k};
            WPlot(seqWeight_1,1);
            TITLE1 = sprintf('%d%s%s',p.L_time,'s-Seq: ',seqName);
            TITLE2 = sprintf('Final activated cells, n=: %d',length(prcInd_confirm));
            title({TITLE1;TITLE2},'Interpreter','None')

            eventStartTs = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{blockInd(j)}).eventStartTs{m};
            eventEndTs = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{blockInd(j)}).eventEndTs{m};
            eventHPeak = seqNMF_DelayDetect_Event_Extract.(seqName).(blockName{blockInd(j)}).eventHPeak{m};


            % delay events
            subplot(4,4,4*(k-1)+[2:4])
            % delay
            for nn = 1:length(prcInd_confirm)
                cellInd = prcInd_confirm(nn);
                tsp = pixelCluster_imec1.clusterInfo{cellInd}.spkTs;
                tspInd = tsp >= delay_Tstart(m) & tsp <= stemStart(m);
                tsp_plot = tsp(tspInd) - delay_Tstart(m) ;
                spkx_delay = tsp_plot * delaySpeed(m);

                xPoints = [spkx_delay';spkx_delay'];
                ypos = nn;
                yPoints = [ypos+zeros(size(tsp_plot'))-0.3;ypos+zeros(size(tsp_plot'))+0.3];
                if ~isempty(tsp_plot)
                    plot(xPoints,yPoints,'k')
                end
                hold on
            end

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
                            plot(xPoints,yPoints,'Color',colorVal(k,:))
                        end
                        hold on
                    end
                end
            end
            xline(maxT*delaySpeed(m),'r--');
            ylim([0 length(prcInd_confirm)+1])
            set(gca, 'YTick', [1,length(prcInd_confirm)], 'YTickLabel', num2cell([1,length(prcInd_confirm)]), 'TickLength', [0, 0]);
            set(gca,'YDir','reverse');
            XLABEL = sprintf('Delay total (cm), %3.1f s', stemStart(m)-delay_Tstart(m));
            xlabel(XLABEL);
            TITLE1 = 'Delay';
            title(TITLE1);
            axis tight
            xlim([0 600])

        end
    end
end

fprintf('Finished figure 7\n');
