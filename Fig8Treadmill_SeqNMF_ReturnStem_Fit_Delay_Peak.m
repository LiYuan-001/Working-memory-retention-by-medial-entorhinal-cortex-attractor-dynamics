% it needs to be stable in each session
% and stable in 2 session combined
%
function Fig8Treadmill_SeqNMF_ReturnStem_Fit_Delay_Peak(inFile,AnalyzeSes)

close all
% Turn off figure pop-ups only for this code block
set(0, 'DefaultFigureVisible', 'on');

addpath(genpath('C:\Users\Li\Documents\MATLAB\seqNMF-master'));

p.colorNum = 5;
colors = get(groot, 'defaultAxesColorOrder');
colorVal = colors(1:p.colorNum, :);

p.timeBin = 200/10^3;
p.L_time = 5; % unit in second
p.seqThresPct = 99;
p.consBinThres = 6;
p.w_min = 1;
p.min_pyrNum = 10; % here make it low to include more sequences, will increase to 15 for quant
p.pValThres = 0.05;
p.minCorrCell = 9; % because matalb corr doesnt calculate p if < 9

% need to be adjusted to automatic detection later
p.linROI = [33, 46];
p.linChoice = [55, 64];
p.xTickLin = [16 40 50 60 68];
p.xTickLabel = {'Return','Dl','St','Ch','Rw'};

p.distanceBin_return = [50,80,110,140];
p.distanceBin_stem = [50,80];
p.distIncre = 10;
% p.returnBase = 200;
% p.stemChoice = 90;

p.savePlot = 1;
p.writeToFile = 0;

% Read in input information
sessInfo = SessInfoImport_Dual(inFile);

for i = AnalyzeSes(1:end)
    close all
    map_All_imec1  = [];
    map_On_imec1 = [];
    map_Off_imec1  = [];
    map_L_Return_imec1 = [];
    map_R_Return_imec1  = [];
    map_L_Choice_imec1 = [];
    map_R_Choice_imec1  = [];
    
    % load synchrinized pos t
    load(fullfile(sessInfo(i).NIDQ,'Processed', 'posT_Sync_NIDQ.mat'));
    % load cell int pyr label in MEC
    cellTypeFile = fullfile(sessInfo(i).NIDQ,'Processed', 'mec_clusterType.mat');
    load(cellTypeFile);
    % load single pike file
    imec1_file = fullfile(sessInfo(i).NIDQ,'Processed', 'pixelCluster_imec1.mat');
    load(imec1_file);
    % % load arm rate file
    % armRateFile = fullfile(sessInfo(i).NIDQ,'Processed', 'CellRate.mat');
    % load(armRateFile);
    %     % periodic file
    %     periodicFile = fullfile(sessInfo(i).NIDQ,'Processed', 'DelayFire_Distance_AutoCorr_Shuffle.mat');
    %     load(periodicFile);
    % load speed
    delaySpeedFile = fullfile(sessInfo(i).NIDQ,'Processed','Delay_Speed.mat');
    load(delaySpeedFile);
    
    % load sequences
    fileName = sprintf('%s%d%s%d%s','Fig8DelaySeq-',ceil(p.L_time*10^3),'ms-',p.timeBin*10^3,'ms.mat');
    seq_file = fullfile(sessInfo(i).NIDQ,'Processed', fileName);
    load(seq_file);
    
    % load shuffled sequences
    fileName = sprintf('%s%d%s%d%s','Fig8CenterSeq-Shuffle-',ceil(p.L_time*10^3),'ms-',p.timeBin*10^3,'ms.mat');
    seq_file = fullfile(sessInfo(i).NIDQ,'Processed', fileName);
    load(seq_file);
    
    % load extracted seq
    fileName = sprintf('%s%d%s%d%s','seqNMF_DelayDetect_Event_Extract_peak-',ceil(p.L_time*10^3),'ms-',p.timeBin*10^3,'ms.mat');
    seq_file = fullfile(sessInfo(i).NIDQ,'Processed', fileName);
    load(seq_file);
    
    % load delayFire map
    delayFile = fullfile(sessInfo(i).NIDQ,'Processed','Fig8DelayTimeMap_2Session.mat');
    load(delayFile);
    
    % get valid cell ind
    %     pyrLabel_imec1 = mec_clusterType.labelInd;
    pyrLabel_imec1 = Fig8_seq_Delay_detect.imec1.pyrInd;
    avgRate_imec1 = mec_clusterType.avgRate;
    depth_imec1 = mec_clusterType.clusterDepth;
    cellLayer_imec1 = mec_clusterType.layerInd;
    depthLayer_imec1 = mec_clusterType.posInd;
    clusterID_imec1 = mec_clusterType.clusterID;
    clusterNum = mec_clusterType.clusterNum;
    pyr = (pyrLabel_imec1 == 1)';
    
    
    % get average lin map for each block
    sessDirs = sessInfo(i).sessDirs;
    spikeMap = cell(clusterNum,1);
    occupMap = cell(clusterNum,1);
    
    spikeMap_On = cell(clusterNum,1);
    spikeMap_Off = cell(clusterNum,1);
    spikeMap_L_Return = cell(clusterNum,1);
    spikeMap_R_Return = cell(clusterNum,1);
    spikeMap_L_Choice = cell(clusterNum,1);
    spikeMap_R_Choice = cell(clusterNum,1);
    
    occupMap_On = cell(clusterNum,1);
    occupMap_Off = cell(clusterNum,1);
    occupMap_L_Return = cell(clusterNum,1);
    occupMap_R_Return = cell(clusterNum,1);
    occupMap_L_Choice = cell(clusterNum,1);
    occupMap_R_Choice = cell(clusterNum,1);
    
    for j = 1:length(sessDirs)
        % load analyzed map per trial
        load(fullfile(sessInfo(i).NIDQ,'Processed','Position',sessDirs{j},'ratesByECLR_imec1.mat'));
        trialNum = length(ratesByECLR_imec1.ECLR);
        
        % SIMPLE WAY: put all trials together without considering
        % differences
        % Accurate way: match left to left and right to right (return,
        % base, choice, reward)
        for n = 1:trialNum
            
            %             maze_time_Temp = sum(ratesByECLR.occupmap{1,1}([1:32,47:end]));
            %             delay_time_Temp = sum(ratesByECLR.occupmap{1,1}([33:46]));
            %             trial_time_Temp = trial_time_Temp + sum((ratesByECLR.occupmap{1,1}));
            
            for k = 1:clusterNum
                spikeMap{k} = [spikeMap{k};ratesByECLR_imec1.spikemap{k,n}];
                occupMap{k} = [occupMap{k};ratesByECLR_imec1.occupmap{k,n}];
                
                if contains(sessDirs{j},'on')
                    spikeMap_On{k} = [spikeMap_On{k};ratesByECLR_imec1.spikemap{k,n}];
                    occupMap_On{k} = [occupMap_On{k};ratesByECLR_imec1.occupmap{k,n}];
                elseif contains(sessDirs{j},'off')
                    spikeMap_Off{k} = [spikeMap_Off{k};ratesByECLR_imec1.spikemap{k,n}];
                    occupMap_Off{k} = [occupMap_Off{k};ratesByECLR_imec1.occupmap{k,n}];
                else
                    error('Block type wrong')
                end
                
                if ratesByECLR_imec1.ECLR(n) == 1 || ratesByECLR_imec1.ECLR(n) == 4
                    spikeMap_R_Return{k} = [spikeMap_R_Return{k};ratesByECLR_imec1.spikemap{k,n}];
                    occupMap_R_Return{k} = [occupMap_R_Return{k};ratesByECLR_imec1.occupmap{k,n}];
                else
                    spikeMap_L_Return{k} = [spikeMap_L_Return{k};ratesByECLR_imec1.spikemap{k,n}];
                    occupMap_L_Return{k} = [occupMap_L_Return{k};ratesByECLR_imec1.occupmap{k,n}];
                end
                
                if ratesByECLR_imec1.ECLR(n) == 1 || ratesByECLR_imec1.ECLR(n) == 2
                    spikeMap_L_Choice{k} = [spikeMap_L_Choice{k};ratesByECLR_imec1.spikemap{k,n}];
                    occupMap_L_Choice{k} = [occupMap_L_Choice{k};ratesByECLR_imec1.occupmap{k,n}];
                else
                    spikeMap_R_Choice{k} = [spikeMap_R_Choice{k};ratesByECLR_imec1.spikemap{k,n}];
                    occupMap_R_Choice{k} = [occupMap_R_Choice{k};ratesByECLR_imec1.occupmap{k,n}];
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
        
        % treadmill on
        spikeMap_AllTrial = sum(spikeMap_On{k},1);
        occupMap_AllTrial = sum(occupMap_On{k},1);
        map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
        map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
        map_On_imec1 = [map_On_imec1;map_Temp_2];
        
        % treadmill off
        spikeMap_AllTrial = sum(spikeMap_Off{k},1);
        occupMap_AllTrial = sum(occupMap_Off{k},1);
        map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
        map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
        map_Off_imec1 = [map_Off_imec1;map_Temp_2];
        
%         % left return
%         spikeMap_AllTrial = sum(spikeMap_L_Return{k},1);
%         occupMap_AllTrial = sum(occupMap_L_Return{k},1);
%         map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
%         map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
%         map_L_Return_imec1 = [map_L_Return_imec1;map_Temp_2];
%         
%         % right return
%         spikeMap_AllTrial = sum(spikeMap_R_Return{k},1);
%         occupMap_AllTrial = sum(occupMap_R_Return{k},1);
%         map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
%         map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
%         map_R_Return_imec1 = [map_R_Return_imec1;map_Temp_2];
%         
%         
%         % left choice
%         spikeMap_AllTrial = sum(spikeMap_L_Choice{k},1);
%         occupMap_AllTrial = sum(occupMap_L_Choice{k},1);
%         map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
%         map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
%         map_L_Choice_imec1 = [map_L_Choice_imec1;map_Temp_2];
%         
%         % right choice
%         spikeMap_AllTrial = sum(spikeMap_R_Choice{k},1);
%         occupMap_AllTrial = sum(occupMap_R_Choice{k},1);
%         map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
%         map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
%         map_R_Choice_imec1 = [map_R_Choice_imec1;map_Temp_2];
        
    end
    
    % get each phase names (no delay etc)
    sessDirs2 = {'on10','off10','on30','off30'};
    for j = 1:length(sessDirs2)
        
        timeMap_Def1.(sessDirs2{j}) = [];
        timeMap_Def2.(sessDirs2{j}) = [];
        pre_Map_Def1.(sessDirs2{j}) = [];
        pre_Map_Def2.(sessDirs2{j}) = [];
        post_Map_Def1.(sessDirs2{j}) = [];
        post_Map_Def2.(sessDirs2{j}) = [];
        
        for k = 1:clusterNum
            %             ind = idx(k);
            ind = k;
            timeMap_Def1.(sessDirs2{j}) = [timeMap_Def1.(sessDirs2{j});Fig8DelayTimeMap_2Session.imec1.(sessDirs2{j}).spikeRate1_Combined_Smooth{ind}];
            timeMap_Def2.(sessDirs2{j}) = [timeMap_Def2.(sessDirs2{j});Fig8DelayTimeMap_2Session.imec1.(sessDirs2{j}).spikeRate2_Combined_Smooth{ind}];
            
            pre_Map_Def1.(sessDirs2{j}) = [pre_Map_Def1.(sessDirs2{j});Fig8DelayTimeMap_2Session.imec1.(sessDirs2{j}).pre_spikeRate1_Combined_Smooth{ind}];
            pre_Map_Def2.(sessDirs2{j}) = [pre_Map_Def2.(sessDirs2{j});Fig8DelayTimeMap_2Session.imec1.(sessDirs2{j}).pre_spikeRate2_Combined_Smooth{ind}];
            
            post_Map_Def1.(sessDirs2{j}) = [post_Map_Def1.(sessDirs2{j});Fig8DelayTimeMap_2Session.imec1.(sessDirs2{j}).post_mapMean_Def1_2Session{ind}];
            post_Map_Def2.(sessDirs2{j}) = [post_Map_Def2.(sessDirs2{j});Fig8DelayTimeMap_2Session.imec1.(sessDirs2{j}).post_mapMean_Def2_2Session{ind}];
        end
    end
    
    % get weight dist for all pyr cells
    shuffle_pyr = Fig8_seq_Center_detect_Shuffle.imec1.pyrInd;
    if any(pyrLabel_imec1 ~= shuffle_pyr)
        error('Error in cell IDs')
    end
    
    seqNum_Shuffle = Fig8_seq_Center_detect_Shuffle.K;
    shuffleTimes = Fig8_seq_Center_detect_Shuffle.shuffleTimes;
    shuffleWeight = zeros(sum(pyr),shuffleTimes*seqNum_Shuffle);
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
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).ratID,'-day',sessInfo(i).Date,'\SeqNMF_ReturnStemFit_DelayPeak-',p.seqThresPct,'-',p.L_time,'s');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    
    
    % get average lin map for each block
    sessDirs = sessInfo(i).sessDirs;
    
    % start taking cells from sequences
    seqNum = Fig8_seq_Delay_detect.K;
    seqWght_pyr = Fig8_seq_Delay_detect.imec1.zpyr.W;
    seqStrength = Fig8_seq_Delay_detect.imec1.zpyr.H_alltime;
    binTime = Fig8_seq_Delay_detect.timeBin_all;
    indSort_pyr = Fig8_seq_Delay_detect.imec1.zpyr.indSort;
    seqWght_pyr_sort = seqWght_pyr(indSort_pyr,:,:);
    
    
    seqNMF_DelayPeak_ReturnStemFit.rat = sessInfo(i).ratID;
    seqNMF_DelayPeak_ReturnStemFit.day = sessInfo(i).Date;
    seqNMF_DelayPeak_ReturnStemFit.timeBin = p.timeBin;
    seqNMF_DelayPeak_ReturnStemFit.L_time = p.L_time;
    seqNMF_DelayPeak_ReturnStemFit.seqThresPct = p.seqThresPct;
    seqNMF_DelayPeak_ReturnStemFit.distIncre = p.distIncre;
    %     seqNMF_ReturnStemFit.returnBase = p.returnBase;
    %     seqNMF_ReturnStemFit.stemChoice = p.stemChoice;
    
    
    for k = 4:seqNum
        seqName = sprintf('seq%d',k);
        
        [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(seqWght_pyr(:,k,:),1);
        indSort_2 = hybrid(:,3);
        seqWght_pyr_sort_2 = seqWght_pyr(indSort_2,:,:);
        w_thres_sort = w_thres2(indSort_2);
        
        selectSeqW = squeeze(seqWght_pyr_sort_2(:,k,:));
        %         selectSeqW(:,[1,2]) = 0;
        %         selectSeqW(:,end-1:end) = 0;
        [selectSeqW_sum, selectSeqW_pos] = max(selectSeqW,[],2);
        pyrInd = find(pyrLabel_imec1 == 1);
        pyrID = clusterID_imec1(pyrInd);
        depth_Pyr = depth_imec1(pyrInd); % pyr cell depth in original order
        pyrID_sort = pyrID(indSort_2); % sorted cluster ID
        pyrIndSort = pyrInd(indSort_2); % sorted cell ind in all cells
        depth_Pyr_sort = depth_Pyr(indSort_2); % sorted pyr depth
        
        pyr_sort_Select_temp = find(selectSeqW_sum > w_thres_sort); % cell ind in sorted way
        % remove cells's position were the same at the sorting
        peakPos = selectSeqW_pos(pyr_sort_Select_temp);
        remainInd = removeConsBin(peakPos,p.consBinThres);
        
        PyrIdx_sort_select_temp = pyrID_sort(pyr_sort_Select_temp); % cluster ID
        pyrInd_sort_select_temp = pyrIndSort(pyr_sort_Select_temp); % cell Ind in all cells
        depth_Pyr_sort_select_temp = depth_Pyr_sort(pyr_sort_Select_temp);
        pyrSelect_temp = indSort_2(pyr_sort_Select_temp);
        selectSeqW_final_temp = selectSeqW_sum(pyr_sort_Select_temp);
        
        seq_cellNumTemp = length(pyr_sort_Select_temp);
        pyrInd_confirm_temp = zeros(seq_cellNumTemp,1);
        for m = 1:seq_cellNumTemp
            pyrInd_confirm_temp(m) = find(clusterID_imec1 == PyrIdx_sort_select_temp(m));
        end
        if any(pyrInd_confirm_temp ~= pyrInd_sort_select_temp)
            error('Sorting is wrong')
        end
        
        pyr_sort_Select = pyr_sort_Select_temp(remainInd);
        PyrId_sort_select = PyrIdx_sort_select_temp(remainInd);
        pyrInd_sort_select = pyrInd_sort_select_temp(remainInd);
        pyrSelect = pyrSelect_temp(remainInd);
        pyrInd_confirm = pyrInd_confirm_temp(remainInd);
        
        % confirmation
        seq_cellNum = length(pyrInd_confirm);
        
        h = figure(2);
        h.Position = [100,100,900,900];
        subplot(2,1,1)
        WPlot(seqWght_pyr_sort(:,:,:),k);
        TITLE1 = sprintf('%s%d%s%s%s%d%s%d','Rat-',sessInfo(i).ratID,'-Day',sessInfo(i).Date,'-',p.L_time,'s-Seq-Delay-Pyr: ',k);
        TITLE2 = 'Sort by all seqs';
        title({TITLE1;TITLE2},'Interpreter','None')
        subplot(2,4,5)
        WPlot(seqWght_pyr_sort_2(:,k,:),1);
        TITLE1 = sprintf('%d%s%d',p.L_time,'s-Seq: ',k);
        TITLE2 = 'Sort by all pyr cells';
        title({TITLE1;TITLE2},'Interpreter','None')
        subplot(2,4,6)
        plot(selectSeqW_sum,length(selectSeqW_sum):-1:1)
        hold on;
        plot(w_thres_sort,length(selectSeqW_sum):-1:1)
        
        subplot(2,4,7)
        WPlot(seqWght_pyr_sort_2(pyr_sort_Select_temp,k,:),1);
        hold on
        plot(peakPos(~remainInd),length(remainInd) - find(remainInd==0),'rx')
        TITLE1 = sprintf('%s%d%s%s%s','Seq: ',k);
        TITLE2 = sprintf('Activated cells, n=: %d',length(remainInd));
        title({TITLE1;TITLE2},'Interpreter','None')
        
        subplot(2,4,8)
        WPlot(seqWght_pyr_sort_2(pyr_sort_Select,k,:),1);
        TITLE1 = sprintf('%s%d%s%s%s','Seq: ',k);
        TITLE2 = sprintf('Final activated cells, n=: %d',length(pyr_sort_Select));
        title({TITLE1;TITLE2},'Interpreter','None')
        
        figName = sprintf('%s%s%d%s%s%s%d%s%d%s',savedir,'\Rat-',sessInfo(i).ratID,'-Day-',sessInfo(i).Date,'-',p.L_time,'s-Seq-',k,'-SeqCellSort');
        print(figName,'-dpng','-r300');
        print(figName,'-depsc', '-vector');
        
        % plot by selected cells in selected sequences
        % indSort_pyr
        h = figure(3);
        h.Position = [100,100,1600,900];
        subplot(2,5,1)
        WPlot(seqWght_pyr_sort_2(pyr_sort_Select,k,:),1);
        yl = ylim;
        TITLE1 = sprintf('%s%d%s%s%s','Rat-',sessInfo(i).ratID,'-Day',sessInfo(i).Date,'Seq-Delay-Pyr');
        TITLE2 = sprintf('%d%s%d%s',p.L_time,'s-Seq-',k,'Selected Pyr Cells');
        title({TITLE1;TITLE2},'Interpreter','None')
        
        for j = 1:length(sessDirs2)
            subplot(2,5,j+1)
            timeMap_Pyr = timeMap_Def1.(sessDirs2{j});
            pre_Map_Pyr = pre_Map_Def1.(sessDirs2{j});
            post_Map_Pyr = post_Map_Def1.(sessDirs2{j});
            cellMapTemp = timeMap_Pyr(pyrInd_confirm,:);
            pre_map1_sort = pre_Map_Pyr(pyrInd_confirm,:);
            post_map1_sort = post_Map_Pyr(pyrInd_confirm,:);
            cellMapTemp_Sort_Norm = [pre_map1_sort,cellMapTemp,post_map1_sort]./max(cellMapTemp,[],2);
            imagesc(flipud(cellMapTemp_Sort_Norm))
            colormap(jet)
            title(sessDirs2{j},'Interpreter','None')
            axis on
            set(gca, 'ydir', 'normal')
            set(gca, 'xtick', [1,size(pre_map1_sort,2),size([pre_map1_sort,cellMapTemp],2),size(cellMapTemp_Sort_Norm,2)]);
            if contains(sessDirs2{j},'10')
                set(gca, 'xticklabels', [-3 0 10 13]);
            else
                set(gca, 'xticklabels', [-3 0 30 33]);
            end
            clim([0 1])
            ylim(yl);
            
        end
        
        subplot(2,5,6)
        linMapTemp = map_All_imec1;
        linMapTemp_Sort = linMapTemp(pyrInd_confirm,:);
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
        
%         subplot(2,5,7)
%         linMapTemp = map_L_Return_imec1;
%         linMapTemp_Sort = linMapTemp(pyrInd_confirm,:);
%         linMapTemp_Sort_Norm_on = linMapTemp_Sort./max(linMapTemp_Sort,[],2);
%         imagesc(flipud(linMapTemp_Sort_Norm_on))
%         set(gca, 'ydir', 'normal')
%         colormap(jet)
%         set(gca, 'xtick', p.xTickLin);
%         set(gca, 'xticklabels', p.xTickLabel);
%         hold on; vertplot(p.linROI, 0, size(linMapTemp_Sort,1), 'r--');
%         vertplot(p.linChoice, 0, size(linMapTemp_Sort,1), 'r--');
%         TITLE1 = sprintf('%s','L_return');
%         title({TITLE1},'Interpreter','None')
%         ylim(yl);
%         
%         subplot(2,5,8)
%         linMapTemp = map_R_Return_imec1;
%         linMapTemp_Sort = linMapTemp(pyrInd_confirm,:);
%         linMapTemp_Sort_Norm_on = linMapTemp_Sort./max(linMapTemp_Sort,[],2);
%         imagesc(flipud(linMapTemp_Sort_Norm_on))
%         set(gca, 'ydir', 'normal')
%         colormap(jet)
%         set(gca, 'xtick', p.xTickLin);
%         set(gca, 'xticklabels', p.xTickLabel);
%         hold on; vertplot(p.linROI, 0, size(linMapTemp_Sort,1), 'r--');
%         vertplot(p.linChoice, 0, size(linMapTemp_Sort,1), 'r--');
%         TITLE1 = sprintf('%s','R_return');
%         title({TITLE1},'Interpreter','None')
%         ylim(yl);
%         
%         subplot(2,5,9)
%         linMapTemp = map_L_Choice_imec1;
%         linMapTemp_Sort = linMapTemp(pyrInd_confirm,:);
%         linMapTemp_Sort_Norm_on = linMapTemp_Sort./max(linMapTemp_Sort,[],2);
%         imagesc(flipud(linMapTemp_Sort_Norm_on))
%         set(gca, 'ydir', 'normal')
%         colormap(jet)
%         set(gca, 'xtick', p.xTickLin);
%         set(gca, 'xticklabels', p.xTickLabel);
%         hold on; vertplot(p.linROI, 0, size(linMapTemp_Sort,1), 'r--');
%         vertplot(p.linChoice, 0, size(linMapTemp_Sort,1), 'r--');
%         TITLE1 = sprintf('%s','L_choice');
%         title({TITLE1},'Interpreter','None')
%         ylim(yl);
%         
%         subplot(2,5,10)
%         linMapTemp = map_R_Choice_imec1;
%         linMapTemp_Sort = linMapTemp(pyrInd_confirm,:);
%         linMapTemp_Sort_Norm_on = linMapTemp_Sort./max(linMapTemp_Sort,[],2);
%         imagesc(flipud(linMapTemp_Sort_Norm_on))
%         set(gca, 'ydir', 'normal')
%         colormap(jet)
%         set(gca, 'xtick', p.xTickLin);
%         set(gca, 'xticklabels', p.xTickLabel);
%         hold on; vertplot(p.linROI, 0, size(linMapTemp_Sort,1), 'r--');
%         vertplot(p.linChoice, 0, size(linMapTemp_Sort,1), 'r--');
%         TITLE1 = sprintf('%s','R_choice');
%         title({TITLE1},'Interpreter','None')
%         ylim(yl);
        
        figName = sprintf('%s%s%d%s%s%s%d%s%d%s',savedir,'\Rat-',sessInfo(i).ratID,'-Day-',sessInfo(i).Date,'-',p.L_time,'s-Seq-',k,'-Selected-Pyr-Temporal&LinMap');
        print(figName,'-dpng','-r300');
        print(figName,'-depsc', '-painters');
        
        % plot by all cells selected sequences
        % indSort_pyr
        h = figure(4);
        h.Position = [100,100,1600,900];
        subplot(2,5,1)
        WPlot(seqWght_pyr_sort_2(:,k,:),1);
        yl = ylim;
        TITLE1 = sprintf('%s%d%s%s%s','Rat-',sessInfo(i).ratID,'-Day',sessInfo(i).Date,'Seq-Delay-Pyr');
        TITLE2 = sprintf('%d%s%d%s',p.L_time,'s-Seq-',k,'All Pyr Cells');
        title({TITLE1;TITLE2},'Interpreter','None')
        
        for j = 1:length(sessDirs2)
            
            subplot(2,5,j+1)
            timeMap_Pyr = timeMap_Def1.(sessDirs2{j});
            pre_Map_Pyr = pre_Map_Def1.(sessDirs2{j});
            post_Map_Pyr = post_Map_Def1.(sessDirs2{j});
            cellMapTemp = timeMap_Pyr(pyrIndSort,:);
            pre_map1_sort = pre_Map_Pyr(pyrIndSort,:);
            post_map1_sort = post_Map_Pyr(pyrIndSort,:);
            cellMapTemp_Sort_Norm = [pre_map1_sort,cellMapTemp,post_map1_sort]./max(cellMapTemp,[],2);
            
            imagesc(flipud(cellMapTemp_Sort_Norm))
            colormap(jet)
            title(sessDirs2{j},'Interpreter','None')
            axis on
            set(gca, 'ydir', 'normal')
            set(gca, 'xtick', [1,size(pre_map1_sort,2),size([pre_map1_sort,cellMapTemp],2),size(cellMapTemp_Sort_Norm,2)]);
            if contains(sessDirs2{j},'10')
                set(gca, 'xticklabels', [-3 0 10 13]);
            else
                set(gca, 'xticklabels', [-3 0 30 33]);
            end
            clim([0 1])
            ylim(yl);
        end
        
        subplot(2,5,6)
        linMapTemp = map_All_imec1;
        linMapTemp_Sort = linMapTemp(pyrIndSort,:);
        linMapTemp_Sort_Norm_on = linMapTemp_Sort./max(linMapTemp_Sort,[],2);
        imagesc(flipud(linMapTemp_Sort_Norm_on))
        colormap(jet)
        set(gca, 'ydir', 'normal')
        set(gca, 'xtick', p.xTickLin);
        set(gca, 'xticklabels', p.xTickLabel);
        hold on; vertplot(p.linROI, 0, size(linMapTemp_Sort,1), 'r--');
        vertplot(p.linChoice, 0, size(linMapTemp_Sort,1), 'r--');
        TITLE1 = sprintf('%s','LinMap_All');
        title({TITLE1},'Interpreter','None')
        ylim(yl);
        
%         subplot(2,5,7)
%         linMapTemp = map_L_Return_imec1;
%         linMapTemp_Sort = linMapTemp(pyrIndSort,:);
%         linMapTemp_Sort_Norm_on = linMapTemp_Sort./max(linMapTemp_Sort,[],2);
%         imagesc(flipud(linMapTemp_Sort_Norm_on))
%         colormap(jet)
%         set(gca, 'ydir', 'normal')
%         set(gca, 'xtick', p.xTickLin);
%         set(gca, 'xticklabels', p.xTickLabel);
%         hold on; vertplot(p.linROI, 0, size(linMapTemp_Sort,1), 'r--');
%         vertplot(p.linChoice, 0, size(linMapTemp_Sort,1), 'r--');
%         TITLE1 = sprintf('%s','L_return');
%         title({TITLE1},'Interpreter','None')
%         ylim(yl);
%         
%         subplot(2,5,8)
%         linMapTemp = map_R_Return_imec1;
%         linMapTemp_Sort = linMapTemp(pyrIndSort,:);
%         linMapTemp_Sort_Norm_on = linMapTemp_Sort./max(linMapTemp_Sort,[],2);
%         imagesc(flipud(linMapTemp_Sort_Norm_on))
%         colormap(jet)
%         set(gca, 'ydir', 'normal')
%         set(gca, 'xtick', p.xTickLin);
%         set(gca, 'xticklabels', p.xTickLabel);
%         hold on; vertplot(p.linROI, 0, size(linMapTemp_Sort,1), 'r--');
%         vertplot(p.linChoice, 0, size(linMapTemp_Sort,1), 'r--');
%         TITLE1 = sprintf('%s','R_return');
%         title({TITLE1},'Interpreter','None')
%         ylim(yl);
%         
%         subplot(2,5,9)
%         linMapTemp = map_L_Choice_imec1;
%         linMapTemp_Sort = linMapTemp(pyrIndSort,:);
%         linMapTemp_Sort_Norm_on = linMapTemp_Sort./max(linMapTemp_Sort,[],2);
%         imagesc(flipud(linMapTemp_Sort_Norm_on))
%         colormap(jet)
%         set(gca, 'ydir', 'normal')
%         set(gca, 'xtick', p.xTickLin);
%         set(gca, 'xticklabels', p.xTickLabel);
%         hold on; vertplot(p.linROI, 0, size(linMapTemp_Sort,1), 'r--');
%         vertplot(p.linChoice, 0, size(linMapTemp_Sort,1), 'r--');
%         TITLE1 = sprintf('%s','L_choice');
%         title({TITLE1},'Interpreter','None')
%         ylim(yl);
%         
%         subplot(2,5,10)
%         linMapTemp = map_R_Choice_imec1;
%         linMapTemp_Sort = linMapTemp(pyrIndSort,:);
%         linMapTemp_Sort_Norm_on = linMapTemp_Sort./max(linMapTemp_Sort,[],2);
%         imagesc(flipud(linMapTemp_Sort_Norm_on))
%         colormap(jet)
%         set(gca, 'ydir', 'normal')
%         set(gca, 'xtick', p.xTickLin);
%         set(gca, 'xticklabels', p.xTickLabel);
%         hold on; vertplot(p.linROI, 0, size(linMapTemp_Sort,1), 'r--');
%         vertplot(p.linChoice, 0, size(linMapTemp_Sort,1), 'r--');
%         TITLE1 = sprintf('%s','R_choice');
%         title({TITLE1},'Interpreter','None')
%         ylim(yl);
%         
        figName = sprintf('%s%s%d%s%s%s%d%s%d%s',savedir,'\Rat-',sessInfo(i).ratID,'-Day-',sessInfo(i).Date,'-',p.L_time,'s-Seq-',k,'-All-Pyr-Temporal&LinMap');
        print(figName,'-dpng','-r300');
        print(figName,'-depsc', '-painters');
        
        close all
        
        orderVec = 1:length(pyrInd_confirm);
        orderVec_weighted = max(selectSeqW(pyr_sort_Select,:),[],2);
        seqNMF_DelayPeak_ReturnStemFit.(seqName).seqCellNum = seq_cellNum;
        seqNMF_DelayPeak_ReturnStemFit.(seqName).cellInd = pyrInd_confirm;
        seqNMF_DelayPeak_ReturnStemFit.(seqName).cellInd_SortedSeq = pyr_sort_Select;
        seqNMF_DelayPeak_ReturnStemFit.(seqName).pyrSelect = pyrSelect;
        
        if seq_cellNum >= p.min_pyrNum
            if any(pyrInd_confirm ~= seqNMF_DelayDetect_Event_Extract.(seqName).cellInd)
                error('Sorting is wrong')
            end
            
            seqNMF_DelayPeak_ReturnStemFit.(seqName).quant = 1;
            h_thres_1 = seqNMF_DelayDetect_Event_Extract.(seqName).h_thres_1;
            % plot ratster plot for each trial from selected seq NMF cells
            % look for best fit val and position
            for j = [1,3,5,7]
                % plot all cells in the sequences in each trial
                % plot selected cells in the sequence in each trial
                % load analyzed positions
                load(fullfile(sessInfo(i).NIDQ,'Processed','Position',sessDirs{j},'ratesByECLR_imec1.mat'));
                delayZoneFile = fullfile(sessInfo(i).NIDQ,'Processed/Position',sessDirs{j}, 'Fig8DelayZonePos.mat');
                load(delayZoneFile);
%                 stemZoneFile = fullfile(sessInfo(i).NIDQ,'Processed/Position',sessDirs{j}, 'Fig8StemZonePos.mat');
%                 load(stemZoneFile);
                delayFile = fullfile(sessInfo(i).NIDQ,'Processed','Position',sessDirs{j}, 'PathZone.mat');
                load(delayFile);
                pathFile = fullfile(sessInfo(i).NIDQ,'Processed/Position',sessDirs{j}, 'pathData.mat');
                pathData = load(pathFile);
                pathData.t = posT_Sync_NIDQ.(sessDirs{j});
                pathLinFile = fullfile(sessInfo(i).NIDQ,'Processed/Position',sessDirs{j}, 'pathDataLinear.mat');
                pathLinData = load(pathLinFile);
                pathLinData.t = posT_Sync_NIDQ.(sessDirs{j});
                
                
                delay_Tstart = Fig8DelayZonePos.delayPos1.startT_Sync;
                stemStart = Fig8DelayZonePos.delayPos1.endT_Sync; % stem start ~= delay pos end
                returnStart = PathZone.posStartT_Sync.Return;
                baseStart = PathZone.posStartT_Sync.Base;
                baseEnd = PathZone.posEndT_Sync.Base;
%                 choiceStart = Fig8StemZonePos.stemPos.endT; % choice start ~= stem pos end
                choiceEnd = PathZone.posEndT_Sync.Choice;
                rewardEnd = PathZone.posEndT_Sync.Reward;
                trialNum = length(delay_Tstart);
                
                
                if contains(sessDirs{j},'on')
                    delaySpeed = Delay_Speed.(sessDirs{j}).speed;
                else
                    delaySpeed = ones(1,trialNum);
                end
                
                if contains(sessDirs{j},'10')
                    maxT = 10;
                else
                    maxT = 30;
                end
                
                for m = 1:trialNum
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
                    
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).bin_start{m} = bin_start;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).bin_Width{m} = bin_Width;
                    
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
                    TITLE1 = sprintf('%s-%s-Trial: %d%',seqName,sessDirs{j},m);
                    title({TITLE1},'Interpreter','None');
                    
                    % this sequence
                    subplot(2,3,3)
                    WPlot(seqWght_pyr_sort_2(pyr_sort_Select,k,:),1);
                    TITLE1 = sprintf('%d%s%s',p.L_time,'s-Seq: ',seqName);
                    TITLE2 = sprintf('Final activated cells, n=: %d',length(pyr_sort_Select));
                    title({TITLE1;TITLE2},'Interpreter','None')
                    
                    h1 = subplot(2,3,2);
                    seqH = seqNMF_DelayDetect_Event_Extract.(seqName).(sessDirs{j}).seqStrength_Delay{m};
                    seqH_filter = seqNMF_DelayDetect_Event_Extract.(seqName).(sessDirs{j}).seqStrength_Delay_Newfit_filter{m};
                    binTimeTemp = seqNMF_DelayDetect_Event_Extract.(seqName).(sessDirs{j}).filter_binTime{m};
                    
                    eventStartTs = seqNMF_DelayDetect_Event_Extract.(seqName).(sessDirs{j}).eventStartTs{m};
                    eventEndTs = seqNMF_DelayDetect_Event_Extract.(seqName).(sessDirs{j}).eventEndTs{m};
                    eventHPeak = seqNMF_DelayDetect_Event_Extract.(seqName).(sessDirs{j}).eventHPeak{m};
                    
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
                    for nn = 1:length(pyrInd_confirm)
                        cellIndTemp = pyrInd_confirm(nn);
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
                        %                         title('Stem+Choice+Reward')
                    end
                    
                    % delay events
                    subplot(2,3,5)
                    if ~isempty(eventStartTs)
                        for kk = 1:length(eventStartTs)
                            colorInd = rem(kk,p.colorNum)+1;
                            if colorInd == 0
                                colorInd = 3;
                            end
                            for nn = 1:length(pyrInd_confirm)
                                cellInd = pyrInd_confirm(nn);
                                tsp = pixelCluster_imec1.clusterInfo{cellInd}.spkTs;
                                tspInd = tsp >= eventStartTs(kk) & tsp <= eventEndTs(kk);
                                
                                tsp_plot = tsp(tspInd) - delay_Tstart(m) ;
                                spkx_delay = tsp_plot * delaySpeed(m);
                                xPoints = [spkx_delay';spkx_delay'];
                                ypos = nn;
                                yPoints = [ypos+zeros(size(tsp_plot'))-0.3;ypos+zeros(size(tsp_plot'))+0.3];
                                if ~isempty(tsp_plot)
                                    plot(xPoints,yPoints,'Color',colorVal(colorInd,:))
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
                    spike_Count_return = zeros(length(bin_start.time.return),length(pyrInd_confirm));
                    spike_ts_Median_return = zeros(length(bin_start.time.return),length(pyrInd_confirm));
                    spike_ts_Mean_return = zeros(length(bin_start.time.return),length(pyrInd_confirm));
                    spike_ts_MedianOrder_return = zeros(length(bin_start.time.return),length(pyrInd_confirm));
                    spike_ts_MeanOrder_return = zeros(length(bin_start.time.return),length(pyrInd_confirm));
                    spike_ts_Corr_rank_return = nan(length(bin_start.time.return),1);
                    spike_ts_Corr_rank_p_return = nan(length(bin_start.time.return),1);
                    spike_ts_Corr_rank_weighted_return = nan(length(bin_start.time.return),1);
                    
                    spike_dist_Median_return = zeros(length(bin_start.time.return),length(pyrInd_confirm));
                    spike_dist_Mean_return = zeros(length(bin_start.time.return),length(pyrInd_confirm));
                    spike_dist_MedianOrder_return = zeros(length(bin_start.time.return),length(pyrInd_confirm));
                    spike_dist_MeanOrder_return = zeros(length(bin_start.time.return),length(pyrInd_confirm));
                    spike_dist_Corr_rank_return = nan(length(bin_start.time.return),1);
                    spike_dist_Corr_rank_p_return = nan(length(bin_start.time.return),1);
                    spike_dist_Corr_rank_weighted_return = nan(length(bin_start.time.return),1);
                    
                    for mm = 1:length(bin_start.time.return)
                        startT = bin_start.time.return(mm);
                        endT = bin_start.time.return(mm) + bin_Width.time.return(mm);
                        
                        spike_ts_CountTemp = zeros(1,length(pyrInd_confirm));
                        spike_ts_MedianTemp = nan(1,length(pyrInd_confirm));
                        spike_ts_MeanTemp = nan(1,length(pyrInd_confirm));
                        
                        spike_dist_MedianTemp = nan(1,length(pyrInd_confirm));
                        spike_dist_MeanTemp = nan(1,length(pyrInd_confirm));
                        
                        for nn = 1:length(pyrInd_confirm)
                            cellInd = pyrInd_confirm(nn);
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
                        
                        % calculated the correlation
                        % 1. rank order 2. weighted rank order 3. cosine similarity
                        if sum(~isnan(order_ts_median))>= p.minCorrCell
                            
                            % time
                            [spike_ts_Corr_rank_return(mm),spike_ts_Corr_rank_p_return(mm)] = corr(orderVec',order_ts_mean', 'Type', 'Spearman', 'Rows', 'complete');
                            
                            validIdx = ~isnan(order_ts_mean);  % Keep only valid entries
                            v1 = orderVec(validIdx);
                            v2 = order_ts_mean(validIdx);
                            w  = orderVec_weighted(validIdx);  % Use weights only on valid entries
                            
                            diffOrder = v1 - v2;
                            validNum = length(w);  % Number of valid entries
                            
                            if validNum > 1  % Need at least 2 values to compute correlation
                                spike_ts_Corr_rank_weighted_return(mm) = 1 - (6 * sum(w' .* (diffOrder.^2))) / (sum(w) * (validNum^2 - 1));
                            else
                                spike_ts_Corr_rank_weighted_return(mm) = NaN;  % Not enough data
                            end
                            
                            % distance
                            [spike_dist_Corr_rank_return(mm),spike_dist_Corr_rank_p_return(mm)] = corr(orderVec',order_dist_mean', 'Type', 'Spearman', 'Rows', 'complete');
                            
                        end
                    end
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spikeCount{m} = spike_Count_return;
                    % time
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_ts_Median{m} = spike_ts_Median_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_ts_Mean{m} = spike_ts_Mean_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_ts_MedianOrder{m} = spike_ts_MedianOrder_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_ts_MeanOrder{m} = spike_ts_MeanOrder_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_ts_Corr_rank{m} = spike_ts_Corr_rank_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_ts_Corr_rank_p{m} = spike_ts_Corr_rank_p_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_ts_Corr_rank_weighted{m} = spike_ts_Corr_rank_weighted_return;
                    % distance
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_dist_Median{m} = spike_dist_Median_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_dist_Mean{m} = spike_dist_Mean_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_dist_MedianOrder{m} = spike_dist_MedianOrder_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_dist_MeanOrder{m} = spike_dist_MeanOrder_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_dist_Corr_rank{m} = spike_dist_Corr_rank_return;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_dist_Corr_rank_p{m} = spike_dist_Corr_rank_p_return;
                    %                     seqNMF_ReturnStemFit.(seqName).(sessDirs{j}).returnCount.spike_dist_Corr_rank_weighted{m} = spike_dist_Corr_rank_weighted_return;
                    
                    
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
                    for nn = 1:length(pyrInd_confirm)
                        cellIndTemp = pyrInd_confirm(nn);
                        tsp = pixelCluster_imec1.clusterInfo{cellIndTemp}.spkTs;
                        %                         % max corr
                        %                         if ~isempty(maxStartT)
                        %                             tspInd = tsp >= maxStartT & tsp <= maxEndT;
                        %                             tsp_plot = tsp(tspInd);
                        %                             [spkx_return,spky,newTsp,spkPosInd] = spk2posInd(tsp_plot,pathLinData.x,pathLinData.y,pathLinData.t);
                        %                             xPoints = [spkx_return';spkx_return'];
                        %                             ypos = nn;
                        %                             yPoints = [ypos+zeros(size(spkx_return'))-0.3;ypos+zeros(size(spkx_return'))+0.3];
                        %                             if ~isempty(tsp_plot)
                        %                                 plot(xPoints,yPoints,'Color',colorVal(2,:))
                        %                             end
                        %                         end
                        
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
                                plot(xPoints,yPoints,'Color',colorVal(1,:))
                            end
                        end
                    end
                    
                    %                     if ~isempty(maxStartT)
                    %                         TEXT = sprintf('max Corr: %1.2f',maxCorr);
                    %                         text(bin_start.distance.return(maxInd),0,TEXT)
                    %                         TEXT = sprintf('pVal: %.3e',maxCorr_p);
                    %                         text(bin_start.distance.return(maxInd),2,TEXT)
                    %                     end
                    
                    if ~isempty(minStartT)
                        TEXT = sprintf('min Corr: %1.2f',minCorr);
                        text(bin_start.distance.return(minInd),0,TEXT)
                        TEXT = sprintf('pVal: %.3e',minCorr_p);
                        text(bin_start.distance.return(minInd),2,TEXT)
                    end
                    
                    % stem
                    spike_Count_stem = zeros(length(bin_start.time.stem),length(pyrInd_confirm));
                    spike_ts_Median_stem = zeros(length(bin_start.time.stem),length(pyrInd_confirm));
                    spike_ts_Mean_stem = zeros(length(bin_start.time.stem),length(pyrInd_confirm));
                    spike_ts_MedianOrder_stem = zeros(length(bin_start.time.stem),length(pyrInd_confirm));
                    spike_ts_MeanOrder_stem = zeros(length(bin_start.time.stem),length(pyrInd_confirm));
                    spike_ts_Corr_rank_stem = nan(length(bin_start.time.stem),1);
                    spike_ts_Corr_rank_p_stem = nan(length(bin_start.time.stem),1);
                    spike_ts_Corr_rank_weighted_stem = nan(length(bin_start.time.stem),1);
                    
                    spike_dist_Median_stem = zeros(length(bin_start.time.stem),length(pyrInd_confirm));
                    spike_dist_Mean_stem = zeros(length(bin_start.time.stem),length(pyrInd_confirm));
                    spike_dist_MedianOrder_stem = zeros(length(bin_start.time.stem),length(pyrInd_confirm));
                    spike_dist_MeanOrder_stem = zeros(length(bin_start.time.stem),length(pyrInd_confirm));
                    spike_dist_Corr_rank_stem = nan(length(bin_start.time.stem),1);
                    spike_dist_Corr_rank_p_stem = nan(length(bin_start.time.stem),1);
                    spike_dist_Corr_rank_weighted_stem = nan(length(bin_start.time.stem),1);
                    
                    for mm = 1:length(bin_start.time.stem)
                        startT = bin_start.time.stem(mm);
                        endT = bin_start.time.stem(mm) + bin_Width.time.stem(mm);
                        
                        spike_ts_CountTemp = zeros(1,length(pyrInd_confirm));
                        spike_ts_MedianTemp = nan(1,length(pyrInd_confirm));
                        spike_ts_MeanTemp = nan(1,length(pyrInd_confirm));
                        
                        spike_dist_MedianTemp = nan(1,length(pyrInd_confirm));
                        spike_dist_MeanTemp = nan(1,length(pyrInd_confirm));
                        
                        for nn = 1:length(pyrInd_confirm)
                            cellInd = pyrInd_confirm(nn);
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
                        % 1. rank order 2. weighted rank order 3. cosine similarity
                        if sum(~isnan(order_ts_median))>= p.minCorrCell
                            
                            % time
                            [spike_ts_Corr_rank_stem(mm),spike_ts_Corr_rank_p_stem(mm)] = corr(orderVec',order_ts_mean', 'Type', 'Spearman', 'Rows', 'complete');
                            
                            validIdx = ~isnan(order_ts_mean);  % Keep only valid entries
                            v1 = orderVec(validIdx);
                            v2 = order_ts_mean(validIdx);
                            w  = orderVec_weighted(validIdx);  % Use weights only on valid entries
                            
                            diffOrder = v1 - v2;
                            validNum = length(w);  % Number of valid entries
                            
                            if validNum > 1  % Need at least 2 values to compute correlation
                                spike_ts_Corr_rank_weighted_stem(mm) = 1 - (6 * sum(w' .* (diffOrder.^2))) / (sum(w) * (validNum^2 - 1));
                            else
                                spike_ts_Corr_rank_weighted_stem(mm) = NaN;  % Not enough data
                            end
                            
                            % distance
                            [spike_dist_Corr_rank_stem(mm),spike_dist_Corr_rank_p_stem(mm)] = corr(orderVec',order_dist_mean', 'Type', 'Spearman', 'Rows', 'complete');
                            
                        end
                    end
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spikeCount{m} = spike_Count_stem;
                    % time
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_ts_Median{m} = spike_ts_Median_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_ts_Mean{m} = spike_ts_Mean_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_ts_MedianOrder{m} = spike_ts_MedianOrder_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_ts_MeanOrder{m} = spike_ts_MeanOrder_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_ts_Corr_rank{m} = spike_ts_Corr_rank_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_ts_Corr_rank_p{m} = spike_ts_Corr_rank_p_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_ts_Corr_rank_weighted{m} = spike_ts_Corr_rank_weighted_stem;
                    % distance
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_dist_Median{m} = spike_dist_Median_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_dist_Mean{m} = spike_dist_Mean_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_dist_MedianOrder{m} = spike_dist_MedianOrder_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_dist_MeanOrder{m} = spike_dist_MeanOrder_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_dist_Corr_rank{m} = spike_dist_Corr_rank_stem;
                    seqNMF_DelayPeak_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_dist_Corr_rank_p{m} = spike_dist_Corr_rank_p_stem;
                    %                     seqNMF_ReturnStemFit.(seqName).(sessDirs{j}).stemCount.spike_dist_Corr_rank_weighted{m} = spike_dist_Corr_rank_weighted_stem;
                    
                    
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
                    for nn = 1:length(pyrInd_confirm)
                        cellIndTemp = pyrInd_confirm(nn);
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
                                plot(xPoints,yPoints,'Color',colorVal(4,:))
                            end
                        end
                        
                        %                         % min corr
                        %                         if ~isempty(minStartT)
                        %                             tspInd = tsp >= minStartT & tsp <= minEndT;
                        %                             tsp_plot = tsp(tspInd);
                        %                             [spkx_return,spky,newTsp,spkPosInd] = spk2posInd(tsp_plot,pathLinData.x,pathLinData.y,pathLinData.t);
                        %                             xPoints = [spkx_return';spkx_return'];
                        %                             ypos = nn;
                        %                             yPoints = [ypos+zeros(size(spkx_return'))-0.3;ypos+zeros(size(spkx_return'))+0.3];
                        %                             if ~isempty(tsp_plot)
                        %                                 plot(xPoints,yPoints,'Color',colorVal(1,:))
                        %                             end
                        %                         end
                    end
                    
                    if ~isempty(maxStartT)
                        TEXT = sprintf('max Corr: %1.2f',maxCorr);
                        text(bin_start.distance.stem(maxInd),0,TEXT)
                        TEXT = sprintf('pVal: %.3e',maxCorr_p);
                        text(bin_start.distance.stem(maxInd),2,TEXT)
                    end
                    
                    %                     if ~isempty(minStartT)
                    %                         TEXT = sprintf('min Corr: %1.2f',minCorr);
                    %                         text(bin_start.distance.return(minInd),0,TEXT)
                    %                         TEXT = sprintf('pVal: %.3e',minCorr_p);
                    %                         text(bin_start.distance.return(minInd),2,TEXT)
                    %                     end
                    
                    
                    
                    subplot(2,3,4)
                    ylim([0 length(pyrInd_confirm)+1])
                    set(gca, 'YTick', [1,length(pyrInd_confirm)], 'YTickLabel', num2cell([1,length(pyrInd_confirm)]), 'TickLength', [0, 0]);
                    set(gca,'YDir','reverse');
                    XLABEL = sprintf('Return + base (cm), %3.1f s', delay_Tstart(m)-returnStart(m));
                    xlabel(XLABEL);
                    TITLE1 = sprintf('%s-Trial: %d return',sessDirs{j},m);
                    title(TITLE1);
                    axis tight
                    xlim([0 200])
                    
                    subplot(2,3,5)
                    %             xlim([0 190])
                    ylim([0 length(pyrInd_confirm)+1])
                    set(gca, 'YTick', [1,length(pyrInd_confirm)], 'YTickLabel', num2cell([1,length(pyrInd_confirm)]), 'TickLength', [0, 0]);
                    set(gca,'YDir','reverse');
                    XLABEL = sprintf('Delay total (cm), %3.1f s', stemStart(m)-delay_Tstart(m));
                    xlabel(XLABEL);
                    TITLE1 = 'Delay';
                    title(TITLE1);
                    axis tight
                    xlim([0 600])
                    
                    subplot(2,3,6)
                    ylim([0 length(pyrInd_confirm)+1])
                    set(gca, 'YTick', [1,length(pyrInd_confirm)], 'YTickLabel', num2cell([1,length(pyrInd_confirm)]), 'TickLength', [0, 0]);
                    set(gca,'YDir','reverse');
                    XLABEL = sprintf('Stem + choice + reward (cm), %3.1f s', rewardEnd(m)-stemStart(m));
                    xlabel(XLABEL);
                    axis tight
                    xlim([230 430])
                    
                    figName = sprintf('%s%s%d%s%d%s%s%s%d%s%d%s%s%s%d',savedir,'\DelayPeak-SeqNMF-Fit-',p.seqThresPct,'-Rat-',sessInfo(i).ratID,'-Day-',sessInfo(i).Date,'-',p.L_time,'s-Seq',k,'-',sessDirs{j},'Trial-',m);
                    print(figName,'-dpng','-r200');
                    print(figName,'-depsc', '-vector');
                    close all
                end
            end
        else
            seqNMF_DelayPeak_ReturnStemFit.(seqName).quant = 0;
        end
    end
    fileName = sprintf('%s%d%s%d%s','seqNMF_DelayPeak_ReturnStemFit-',ceil(p.L_time*10^3),'ms-',p.timeBin*10^3,'ms.mat');
    save(fullfile(sessInfo(i).NIDQ,'Processed',fileName), 'seqNMF_DelayPeak_ReturnStemFit');
    clear seqNMF_DelayPeak_ReturnStemFit
    fprintf('Finished analysis for session %d\n',i);
end
% Turn on figure pop-ups only for this code block
set(0, 'DefaultFigureVisible', 'on');

end
