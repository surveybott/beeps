%%% From Chris via Gabriela %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gabriela  6:20 PM
% Analyses average trial:
% -      Chris’s comments is that overall the data is way to noisy and the 
% baseline could be better selected
% How to do z scoring and analysis:
% Trim for each participant beginning and end make sure you don’t have 
% spill from before and after
% z scoring whole time series for each participant within experiment
% Step 1: identify for each event a 1 s period before the beep
% Step 2: average that pre-trial period out for each trial and then average
% across trials and across participants
% For plotting
% -      1 sec jitter and the full trial to see the event
% -      Trials are too noisy and the question is why?
% -      So now from this new metric we plot the average response for each 
% experiment and confidence intervals
% For habituation:
% Do a sliding window approach
% first time point has average event (just like before event – immediate 
% baseline before the event) to 4th timepoint average event; then plot 2nd 
% to 5th trial and so on
% 6:23
% For luminance:
% -      Frame by frame changes are way too fast no way that pupil can 
% react but our lux measures too slow
% -      What we can see from our lux measures is that there’s very little 
% variability especially in the beeps, so no correction is needed here
% -      We can do corrections at the end for elf and bubblegum but in 100 
% to 200 ms resolutions after converting RGB values (see paper shared)
% 
% 6:23
% basically we won't care for luminance just yet and address the other 
% analyses first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all;


% path to pspm
addpath('C:/Program Files (x86)/MATLAB/R2013a Student/PsPM/v4.2.1')



%% fetch edf files for external conversion
% edfTag = '*.edf';
% edfSet = dir(edfTag);
% for i = 1:length(edfSet)
%     edfNames{i,1} = edfSet(i,1).name;
% end



%% after test/uncertain file removal, checking asc-mat nfile equivalence
% length(dir('*.asc')) == length(dir('*.mat'))
% ascSet = dir('*.asc');
% matSet = dir('*.mat');
% for i = 1:length(ascSet)
%     str = ascSet(i,1).name;
%     ascNames{i,1} = regexprep(str,'.asc','','ignorecase');
% end
% for i = 1:length(matSet)
%     str = matSet(i,1).name;
%     matNames{i,1} = regexprep(str,'.mat','','ignorecase');
% end
% setdiff(matNames,ascNames)
% setdiff(ascNames,matNames)

%%% Current sample sizes: Beeps (18), Bubblegum (17), Elf (18)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
%              -- Import data via PsPM before continuing. --              %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%



%% pupil preprocessing for all files/conditions
% pspmSet = dir('pspm*');
% options.channel='pupil_l';
% options.channel_action='add';
% options.channel_combine='pupil_r';
% for i = 1:length(pspmSet)
%     pspm_pupil_pp(pspmSet(i,1).name,options);
% end



%% fetch data files - just beeps, adults, control
allFileSet = dir('beeps*AC*.mat');
pspmSet    = dir('*pspm*beeps*AC*.mat'); 



%% setting up variables
averageTrial      = [];
fullBinMeans      = [];
fullPDdist        = [];
sdprops           = [];
avgSegMeans       = [];
avgTrial          = [];
firstTrial        = [];
avgBaseline       = [];
overallHalfOne    = [];
overallHalfTwo    = [];
oneSecBaseline    = [];
avgBaselineFull   = [];
trialDescriptives = [];
meanTrial         = [];
threeTrials       = [];
pupilCorrs        = [];
nanRatios         = [];
allFiveTrials     = [];
numTrials         = [];

fiftyTrials = nan(125000,length(pspmSet));
tenTrials   = nan(25000,length(pspmSet));
fiveTrials  = nan(12500,length(pspmSet));
wholeSeries = nan(400000,length(pspmSet));
firstFive   = nan(2500,length(pspmSet));

CT.header = {'PART','GROUP','COND','BIN','PUPIL'};
CT.descriptives = [];

load('beeps_NT_data.mat')
beeps.adult.sub = {};



%% computing luminance values for reweighting of PD timeseries
load('luminance_values.mat') % relative values from videos
load('lux_data.mat') % luxmeter data

%%% finding framewise proportional change
% diffs = [];
% for i = 1:length(lum.base.beeps)-1
%     diffs(i,1) = lum.base.beeps(i+1)/lum.base.beeps(i);
% end
% diffs(end+1) = diffs(end);



%% visual checking for data/distributions
for i = 1:length(pspmSet)
    
    %%% load files and find demographics
    load(pspmSet(i).name);
    load(allFileSet(i).name);
    Str            = pspmSet(i).name;
    Key1           = '_';
    Key2           = 'A';
    Index1         = strfind(Str, Key1);
    Index2         = strfind(Str, Key2);
    part.id        = sscanf(Str(Index2(1,1)+2:Index2(1,1)+3),'%g',1);
    part.group     = sscanf(Str(Index2(1,1):Index2(1,1)+1),'%s',1);
    part.condition = sscanf(Str(Index1(1,1)+1:Index1(1,1)+3),'%s',1);
    part.nanRatio  = infos.source.chan_stats{3,1}.nan_ratio;
    part.numMarks  = length(data{7,1}.data);
    clear str Key1 Key2 Index1 Index2
    
    %%% plot baseline
    % plot(data{9,1}.data(data{7,1}.data(2)*1000:data{7,1}.data(3)*1000));
    % legend(int2str(part.id))
    
    %%% plot data
    % plot(data{9,1}.data(data{7,1}.data(3)*1000:end));
    % legend(int2str(part.id))
    
    %%% plot distribution
    % hist(data{9,1}.data(data{7,1}.data(2)*1000:data{7,1}.data(3)*1000))
    % legend(int2str(part.id))
    
    % pause()
    
    %%% save image if needed
    % saveas(gcf,['hist_beeps_CC_baseline_' int2str(part.id) '.png'])
    % close all

    %%% collect all participant distributions and mean pupils
    fullPDdist = [fullPDdist;data{9,1}.data(data{7,1}.data(2):data{7,1}.data(3))];
    avgFirstTrial(i) = nanmean(data{9,1}.data(data{7,1}.data(2):data{7,1}.data(3)));
    
end

%%% plot distribution for all participants
% hist(fullPDdist)

%%% compute outlier 3SD proportions and cutoffs
dist.m = nanmean(fullPDdist);
dist.sd = nanstd(fullPDdist);
dist.upProp = sum(fullPDdist>(dist.m+3*dist.sd))/length(fullPDdist);
dist.dnProp = sum(fullPDdist<(dist.m-3*dist.sd))/length(fullPDdist);
dist.tProp = dist.upProp + dist.dnProp;



%% main processing loop
for i = 1:length(pspmSet)
    
    % --- load files and find demographics
    load(pspmSet(i).name);
    load(allFileSet(i).name);
    Str            = pspmSet(i).name;
    Key1           = '_';
    Key2           = 'A';
    Index1         = strfind(Str, Key1);
    Index2         = strfind(Str, Key2);
    part.id        = sscanf(Str(Index2(1,1)+2:Index2(1,1)+3),'%g',1);
    part.group     = sscanf(Str(Index2(1,1):Index2(1,1)+1),'%s',1);
    part.condition = sscanf(Str(Index1(1,1)+1:Index1(1,1)+3),'%s',1);
    clear str Key1 Key2 Index1 Index2

    % --- set up participant metadata
    part.meta.ratios.blink = infos.source.chan_stats{3,1}.blink_ratio;
    part.meta.ratios.saccade = infos.source.chan_stats{3,1}.saccade_ratio;
    part.meta.ratios.nan  = infos.source.chan_stats{3,1}.nan_ratio;
    part.meta.gazeCoords = infos.source.gaze_coords;
    part.meta.samplingRate = data{3,1}.header.sr; % original sampling rate
    part.meta.units = data{3,1}.header.units;
    part.meta.numMarks  = length(data{7,1}.data);
    
    % --- set up pupil data and baseline
    markerSR = int64(data{7,1}.data*1000); % markers x pspm sampling rate
    pupil.baseline.data = data{9,1}.data(1:markerSR(3));
    pupil.baseline.mean = nanmean(pupil.baseline.data);
    pupil.baseline.sd = nanstd(pupil.baseline.data);
    pupil.data = data{9,1}.data;
    pupil.data(1:markerSR(3)) = NaN; % remove basline, retain structure
    pupil.mean = nanmean(pupil.data);
    pupil.sd = nanstd(pupil.data);

    % --- set up gaze data; error - breaks the loop
    % pupil.gazeL(:,1) = data{1,1}.data; % x coords
    % pupil.gazeL(:,2) = data{2,1}.data; % y coords
    % pupil.gazeR(:,1) = data{4,1}.data; % x coords
    % pupil.gazeR(:,2) = data{5,1}.data; % y coords
    
    % --- compute correlations between pupils
    pupilCorrs(i,1) = corr(data{3,1}.data(markerSR(3):end), ...
        data{6,1}.data(markerSR(3):end), 'rows','complete');
    

    
    %% --- trimming outliers
    % pupil.data(pupil.data>(pupil.mean+3*pupil.sd)) = NaN;
    % pupil.data(pupil.data<(pupil.mean-3*pupil.sd)) = NaN;
    
    
    
    %% --- luminance correction attempts
    % ---- lum = 179*(R_irradiance*.265 + G_irradiance*.67 + B_irradiance*.065)
    % ---- https://radiance-online.org/pipermail/radiance-general/2011-July/007954.html
    % ---- to do: extend/interpolate lumvec to same timescale as PD
    
    % lumvec = repmat(lum.trial.beeps,1,floor(length(pupil.data)/4999));
    % lumvec = [lumvec lum.trial.beeps(1:length(pupil)-length(lumvec))];
    % pupil.data = pupil.data./lumvec';
        
    % upsampling baseline luminance vector
    % participants only see video 1.8849 times
    % video is 166s long at 30.1144578313253 f/s
    % need to upsample lum to 1000 f/s, then need length = pupil
    % this isn't perfect, but it's better than original method
    % lumvec = repmat(interp(lum.base.beeps,33),1,2)';
    % lumvec = lumvec(1:length(pupil));
	% pupil.data = pupil.data./lumvec;
    
    
    
    %% --- z-scoring baseline and pupil
    for k = 1:length(pupil.data)
        pupil.data(k) = (pupil.data(k)-pupil.mean)/pupil.sd;
    end
    for k = 1:length(pupil.baseline.data)
        pupil.baseline.data(k) = (pupil.baseline.data(k)- ...
            pupil.baseline.mean)/pupil.baseline.sd;
    end
    
    
    %% --- sliding window attempts (slow)
	% pupilWindow = slidingavg_nan(pupil.data,6000);
	% pupil.data = slidingavg_nan(pupil,6000);
	% pupil.data = sliding_window(pupil,6000);
	% wholeSeries(1:length(pupil.data),i) = sliding_window(pupil.data,100);
    
    
    
    %% --- computing basline and data means
    blines(i,1:2) = [part.id pupil.baseline.mean];
    pMeans(i,1:2) = [part.id pupil.mean];
    
    
    
    %% --- setting up windows for segmenting trial data
    windowStart = [];
    windowEnd   = [];
    for j = 1:length(markerSR)-3
        windowStart(j) = markerSR(j+2);
        windowEnd(j)   = markerSR(j+3);
    end

    
    
    %% --- baseline correction attempts    
    % ---- correcting based on baseline stimulus
    % pupil.corrected = pupil.data-pupil.baseline.mean;
    
    % ---- correcting based on first trial
    % pupil.corrected = pupil.data-nanmean(pupil.data(markerSR(3):markerSR(4)));
    
    % ---- correcting based on mean of last 200ms of previous trial
    % iti2 = [];
    % for k = 1:length(windowEnd)
    %     iti2(:,k) = pupil(windowEnd(k)-200:windowEnd(k));
    % end
    
    
    
    % ---- correcting based on mean of full previous trial; works
    pupil.corrected = nan(length(pupil.data),1);
    pupil.corrected(windowStart(1):windowEnd(1)) = ...
        pupil.data(windowStart(1):windowEnd(1));
    for k = 1:length(windowEnd)-1
        pupil.corrected(windowStart(k+1):windowEnd(k+1)) = ...
            pupil.data(windowStart(k+1):windowEnd(k+1)) - ...
            nanmean(pupil.data(windowStart(k):windowEnd(k)));
    end
    pupil.corrected(windowStart(1):windowEnd(1)) = ...
        pupil.data(windowStart(1):windowEnd(1)) - ...
        nanmean(pupil.data(windowStart(1):windowEnd(1)));

    
    
    % ---- zero-point baseline correction; issues with missing data
    % pupil.corrected = nan(length(pupil.data),1);
    % pupil.corrected(windowStart(1):windowEnd(1)) = ...
    %     pupil.data(windowStart(1):windowEnd(1))-pupil.data(windowStart(1)+1);
    % for k = 1:length(windowEnd)-1
    %     pupil.corrected(windowStart(k+1):windowEnd(k+1)) = ...
    %         pupil.data(windowStart(k+1):windowEnd(k+1)) - ...
    %         pupil.data(windowStart(k)+1);
    % end
    
    

    %% --- segmenting data into trials
    first = {};
    for k = 1:length(windowStart)        
        first{k} = pupil.corrected(windowStart(k):windowEnd(k));        
    end
    % first{end+1}=pupil.corrected(windowEnd(end):end); % tail end of series
    
    
    
    %% --- downsampling segments back to 500Hz
    for k = 1:length(first)
        first{k} = downsample(first{k},2);
    end

    
    
    %% --- comparing baseline with first trial
    baseSegs = []; % segmenting the baseline
    for k = 1:10
        baseSegs(:,k) = pupil.baseline.data(k*1000-1000+1:k*1000-1000+1000);
    end
    
    % ---- compute average 1 second of baseline
    avgBaseSeg = []; 
    for k = 1:1000
        avgBaseSeg(k,1) = nanmean(baseSegs(k,:));
    end
    
    % ---- downsample and gather baseline segments
    avgBaseSeg = downsample(avgBaseSeg,2);
    avgBaseline(:,i) = avgBaseSeg;
    
    % ---- trim and gather first trials 
    firstTrial(:,i) = first{1}(1:500);
    
    % ---- visualizing baseline vs 1st trial for current participant
    % plot(avgBaseSeg,'LineWidth',3)
    % hold on
    % plot(firstTrial(:,i),'LineWidth',3)
    % legend({'Baseline','First Trial'})
    % saveas(gcf,['plot_beeps_CC_baseline_vs_1st_' int2str(partID) '.png'])
    % pause()
    % close all
    
    
    %% --- comparing average first half and average second half trials
    halfOne = [];
    halfTwo = [];
    ind = 1;
    for k = 1:length(first)
        if length(first{k}) >= 500
            halfOne(:,ind) = first{k}(1:500);
            ind = ind + 1;
        end        
    end    
    halfTwo = halfOne(:,floor(length(first)/2):end);    
    halfOne = halfOne(:,1:floor(length(first)/2));
    
    % ---- averaging each half
    avgHalfOne = [];
    for k = 1:length(halfOne(:,1))
        avgHalfOne(k,1) = nanmean(halfOne(k,:));
    end
    avgHalfTwo = [];
    for k = 1:length(halfTwo(:,1))
        avgHalfTwo(k,1) = nanmean(halfTwo(k,:));
    end
    
    % ---- gathering each half across participants
    overallHalfOne(:,i) = avgHalfOne;
    overallHalfTwo(:,i) = avgHalfTwo;
    
    % ---- visualizing two halves for current participant
    % plot(avgHalfOne)
    % hold on
    % plot(avgHalfTwo)
    % legend({'1st Half Trial','2nd Half Trial'})
    % pause()
    % close all


    
    %% --- finding shortest trial, in case needed for trimming
    ls = [];
    for k = 1:length(first)
        ls(k) = length(first{1,k});
    end
    ls = sort(ls);
    
    
    
    %% --- trimming and concatenating segments into matrix for averaging
    trimmed = [];
    for k = 1:length(first)
        trimmed = [trimmed first{1,k}(1:500)];
    end

    
    
    %% --- gathering trials 1 to 5
    firstFive(:,i) = [trimmed(:,1);trimmed(:,2);trimmed(:,3); ...
        trimmed(:,4);trimmed(:,5)];
        
    
    
    %% --- apply sliding window to trial segments
    trimmedFull = [];
    trimmedTen = [];
    for k = 1:length(trimmed(1,:))
        trimmed(:,k) = sliding_window(trimmed(:,k),250);
        trimmedFull = [trimmedFull;trimmed(:,k)];
    end
    
    trimCatch = {};
    if length(first) >= 295
        for k = 1:5
            trimCatch{k} =  trimmed(:,k*59-58:k*59) ;
        end
        
        fiveTrials = [];
        for k = 1:length(trimCatch)
            for j = 1:length(trimCatch{k}(:,1))
                fiveTrials = [fiveTrials;nanmean(trimCatch{k}(j,:))];
            end
        end
        allFiveTrials = [allFiveTrials fiveTrials];
    end
    
    
    %% --- computing average trial
    for k = 1:500
        avgTrial(k,i) = nanmean(trimmed(k,:));
    end
    
    % plot(avgTrial(:,i))
    % xline(55)
    % legend(int2str(part.id),'Beep Offset')
    % ylim([-0.1 0.1])
    % pause()
    % close all
    
    
    
    %% --- gathering into bins
    col = [];
    bins = [];
    for k = 1:length(trimmed(1,:))/6
        col = [];
        for l = 1:6            
            col = [col;trimmed(:,l*k-k+1)];
        end
        bins = [bins col];        
    end
    
    
    
    %% --- averaging
    segMeans = [];
    for l = 1:length(bins(1,:))
        segMeans(l,1) = nanmean(bins(:,l));
    end
    
    % threeTrials = [threeTrials segMeans];
    % for l = 1:floor(length(segMeans)/6)
    %     binMeans(l,1) = nanmean(segMeans(l*6-6+1:l*6-6+6,1));
    % end
    
    
    

    %% --- set outliers to NaN
    % segMeans(segMeans<sgM-sgSD*3|segMeans>sgM+sgSD*3) = NaN;
    % avgSegMeans(:,i) = segMeans;
    
    
    
    %% --- computing and preparing descriptives
    binMeans     = segMeans;
    partNum      = repmat(part.id,length(binMeans),1);
    binNum       = [1:length(binMeans)]';
    groupNum     = repmat(0,length(binMeans),1);
    condNum      = repmat(0,length(binMeans),1);
    fullBinMeans = [fullBinMeans; partNum groupNum condNum binNum binMeans];
    
    
    
    %% --- setting up stripped down struct for Gabriela  
    beeps.adult.sub{i}.zpupil_corrected = pupil.corrected;
    beeps.adult.sub{i}.zbaseline = pupil.baseline;
    beeps.adult.sub{i}.markers = markerSR;
    % saveas()
    
    
    
    %% --- descriptives
    % binMeans = segMeans;
    % partNum = repmat(partID,length(binMeans),1);
    % binNum = [1:length(binMeans)]';
    % groupNum = repmat(0,length(binMeans),1);
    % condNum = repmat(0,length(binMeans),1);
    % fullBinMeans = [fullBinMeans; partNum groupNum condNum binNum binMeans];
    
    % partNum = repmat(partID,900,1);
    % binNum = [1:900]';
    % groupNum = repmat(0,900,1);
    % condNum = repmat(0,900,1);
    % trialDescriptives = [trialDescriptives; partNum groupNum condNum binNum avgTrial(1:900)'];
    
end



avgPupilCorr = nanmean(pupilCorrs);
rngPupilCorr = [nanmin(pupilCorrs) nanmax(pupilCorrs)];



%% --- plotting average trial
% ---- collapse across participants
mResp  = [];
sdResp = [];
for k = 1:length(avgTrial(:,1))
    mResp(k,1)  = nanmean(avgTrial(k,:));
    sdResp(k,1) = nanstd(avgTrial(k,:));
end
mResp   = mResp-mResp(1,1); % set 1st point to zero

% ---- setting up confidence intervals
x       = (1:length(mResp));
Data    = mResp';
Data_sd = sdResp';
x_ax    = x;
X_plot  = [x_ax, fliplr(x_ax)];
Y_plot  = [Data-(1.96.*(Data_sd/sqrt(18))), fliplr(Data+(1.96.*(Data_sd/sqrt(18))))];

% ---- plot data and confidence intervals
hold on 
plot(x_ax, Data, 'blue', 'LineWidth', 3)
fill(X_plot, Y_plot , 1,....
        'facecolor','blue', ...
        'edgecolor','none', ...
        'facealpha', 0.3);
hold off 
set(gca,'FontSize',20)
set(gcf,'position',[10,10,800,800])
hold on
title({'Average Beeps Trial, Adult Control'},'FontSize',30)
xlabel('Time (s)','FontSize',25)
ylabel('z-scored Pupil Size (a.u.)','FontSize',25)
% ylim([-1.6 0.6])
xticks([0 125 250 375 500])
xticklabels({'0','0.25','0.5','0.75','1'})


