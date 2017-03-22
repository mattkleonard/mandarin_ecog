subj = 'EC143';
blocks = {'B21','B22','B23','B24','B26','B27','B32','B33'};
% behav_blocks = [1 2 3 4 5 6 ];
behav_ecog_correspondence = [1 21 ; 2 22 ; 3 23 ; 4 24 ; 5 26 ; 6 27 ; 7 32 ; 8 32 ; 9 32 ; 10 33 ; 11 33 ; 12 33];
hemi = 'lh';

rootdir = '/Users/mattleonard/Documents/Research/pia/data_store1/human/prcsd_data';

% flags
tanh_flag = 0;                  % whether to perform tahn transform to fix outliers
find_events_flag = 1;           % whether to find events in ANIN channels
make_dat_struct_flag = 1;       % whether to make dat structure
    prestim_zscore_flag = 1;        % whether to z-score relative to pre-stim
plot_all_elecs_flag = 0;        % whether to plot average ERPs for all channels
    plot_raster_flag = 0;           % whether to plot single trial rasters

% time variables
timeLim = [1 1];                        % time limits (sec) of epoching
taxis = -timeLim(1):1/100:timeLim(2);   % create time axis
prestim_zscore_time = [-0.5 0];         % time window to use for pre-stim z-score

% plotting variables
block_nums = [1:6];
syls = {'bu','di','ma'}; % {'bu','di','lu','ma','mi'};
spkrs = {'a','i'};
tones = [1 2 3 4]; % [1 2 3 4];
acc = [0 1]; % [0 1];

% colors
clr = [0.667, 0.224, 0.224 ; 0.176, 0.533, 0.176 ; 0.133, 0.40, 0.40 ; 0.667, 0.424, 0.224 ; 0.863 0.502 0.698];

%% REMOVE ARTIFACTS USING TANH

if tanh_flag
    for b = 1:length(blocks)
        dirlist = dir([rootdir '/' subj '/' subj '_' blocks{b} '/HilbAA_70to150_8band_notanh/']);
        dirlist = dirlist(~[dirlist.isdir]);
        
        for i = 1:length(dirlist)
            fprintf('Channel [%d] of [%d]\n',i,length(dirlist));
            [d,fs] = readhtk([rootdir '/' subj '/' subj '_' blocks{b} '/HilbAA_70to150_8band_notanh/' dirlist(i).name]);
            d = mean(d,1);
            mm1 = 10*std2(log10(d));
            mm2 = mm1*tanh(log10(d)/mm1);
            writehtk([rootdir '/' subj '/' subj '_' blocks{b} '/HilbAA_70to150_8band/' dirlist(i).name],mm2,fs);
        end
    end
end

%% FIND EVENTS

if find_events_flag
    
    load([rootdir '/../../../../data/' subj '/mandarin/behavior/' subj '_behav.mat']);
    allEvents = [];
    
    for b = 1:length(blocks)
        dirlist = dir([rootdir '/' subj '/' subj '_' blocks{b} '/HilbAA_70to150_8band/']);
        dirlist = dirlist(~[dirlist.isdir]);
        
        evnts = [];
%         load([rootdir '/' subj '/' subj '_' blocks{b} '/Artifacts/badTimeSegments.mat']);
%         load([rootdir '/' subj '/' subj '_' blocks{b} '/Analog/allEventTimes.mat']);
        
        % get stimuli for this ECoG block
        stimlist = [];
        behav_blocks = find(behav_ecog_correspondence(:,2) == str2double(blocks{b}(2:end)));
        block_idx = find(ismember(behav.block,behav_blocks));
        
        for i = 1:length(block_idx)
            stimlist{i} = [behav.syllable{block_idx(i)} num2str(behav.tone(block_idx(i))) '-' behav.speaker{block_idx(i)} 'N'];
        end
        
        evnt = DetectEventsQuick_phonrest([rootdir '/' subj '/' subj '_' blocks{b}],...
            '/Users/mattleonard/Documents/Research/tasks/ToneCat_ptb3/FindEvents',...
            {[subj '_' blocks{b}]},...
            stimlist,...
            stimlist,...
            2,...
            [],...
            1);
        
        evnts = [evnts evnt];
        allEvents = [allEvents evnts];
        
        save([rootdir '/' subj '/' subj '_' blocks{b} '/Analog/evnts.mat'],'evnts','-v7.3');
    end
end

%% CREATE DAT STRUCTURE

if make_dat_struct_flag
    dat.subj = subj;
    dat.hemi = hemi;
    dat.behav = behav;
    
    for b = 1:length(blocks)
        dirlist = dir([rootdir '/' subj '/' subj '_' blocks{b} '/HilbAA_70to150_8band/']);
        dirlist = dirlist(~[dirlist.isdir]);
        
        load([rootdir '/' subj '/' subj '_' blocks{b} '/Analog/evnts.mat']);
        
        fid = fopen([rootdir '/' subj '/' subj '_' blocks{b} '/Artifacts/badChannels.txt']);
        badChans = textscan(fid,'%d');
        badChans = badChans{1};
        fclose(fid);
        
        badTrials = Find_bad_trials(rootdir,subj,blocks{b},[evnts.StartTime],timeLim);
        
        [~,fs] = readhtk([rootdir '/' subj '/' subj '_' blocks{b} '/HilbAA_70to150_8band/' dirlist(1).name]);
        dat.hg_stim = NaN(length(dirlist),...
            length(round((evnts(1).StartTime * fs/4) - (timeLim(1) * fs/4)):round((evnts(1).StartTime * fs/4) + (timeLim(2) * fs/4))),...
            length(evnts) - length(badTrials));
        
        for i = 1:length(dirlist)
            fprintf('Channel [%d] of [%d]\n',i,length(dirlist));
            [d,fs] = readhtk([rootdir '/' subj '/' subj '_' blocks{b} '/HilbAA_70to150_8band/' dirlist(i).name]);
            d = mean(d,1);
            
            if ~prestim_zscore_flag
                d = (d - mean(d)) / std(d);
            end
            
            d = resample(d,1,4);
            fs = fs / 4;
            
            k = 1;
            for j = 1:size(evnts,2)
                if ~ismember(j,badTrials)
                    dat.hg_stim(i,:,k) = d(round((evnts(j).StartTime * fs) - (timeLim(1) * fs)):round((evnts(j).StartTime * fs) + (timeLim(2) * fs)));
                    if prestim_zscore_flag
                        dat.hg_stim(i,:,k) = (dat.hg_stim(i,:,k)...
                            - nanmean(dat.hg_stim(i,...
                            find(taxis == prestim_zscore_time(1)):find(taxis == prestim_zscore_time(2)),...
                            k),2)) / ...
                            nanstd(dat.hg_stim(i,...
                            find(taxis == prestim_zscore_time(1)):find(taxis == prestim_zscore_time(2)),...
                            k),[],2);
                    end
                    k = k + 1;
                end
            end
        end
    end
    dat.gridOrient = [];
    dat.badChans = badChans;

    save([rootdir '/' subj '/' subj '_dat.mat'],'dat','-v7.3');
end

%%
%% PLOTTING

if plot_all_elecs_flag
    
    if ~exist('dat','var')
        load([rootdir '/../../../userdata/matt_l/mandarin/' subj '/data/' subj '_dat.mat']);
    end
    
    %% GET TRIALS
    
    trls = [];
    for i = 1:length(tones)
        trls{i} = find((ismember(dat.behav.block,block_nums)) & ...
            (ismember(dat.behav.tone,tones(i))));
    end
    %     for i = 1:length(syls)
    %         trls{i} = find((ismember(dat.behav.block,block_nums)) & ...
    %             (strcmpi(dat.behav.syllable,syls{i})));
    %     end
    
    for i = 1:length(trls)
        for j = 1:length(trls{i})
            for k = 1:size(dat.hg_stim,1)
                if ~ismember(i,dat.badChans)
                    if ~isempty(find(dat.hg_stim(k,:,trls{i}(j)) > 20))
                        badTrials(k,i,j) = 1;
                    else
                        badTrials(k,i,j) = 0;
                    end
                end
            end
        end
    end
    
    for i = 1:length(trls)
        for j = length(trls{i}):-1:1
            if any(badTrials(:,i,j),1)
                trls{i}(j) = [];
            end
        end
    end
    %%
    gridOrient = dat.gridOrient';
    
    figure;
    for i = 1:size(dat.hg_stim,1)
        if ~ismember(i,dat.badChans)
            for j = 1:length(trls)
                tmp_ylim(i,:) = [min(min(squeeze(nanmean(dat.hg_stim(i,:,trls{j}),3)))) max(max(squeeze(nanmean(dat.hg_stim(i,:,trls{j}),3))))];
            end
        end
    end
    ylims = [min(min(tmp_ylim)) max(max(tmp_ylim))];
    
    
    for i = 1:size(gridOrient,1)*size(gridOrient,2)
        if ~isnan(gridOrient(i))
            p = plotGridPosition_new(i,size(gridOrient,1)*size(gridOrient,2),ceil(size(gridOrient,1)));
            subplot('Position',p);
            
            for j = 1:length(trls)
                shadedErrorBar(taxis,...
                    squeeze(nanmean(dat.hg_stim(gridOrient(i),:,trls{j}),3)),...
                    squeeze(ste(dat.hg_stim(gridOrient(i),:,trls{j}),3)),...
                    {'-','Color',clr(j,:),'LineWidth',1},1);
                hold on;
            end
            axis tight;
            
            set(gca,'YLim',ylims,'XTickLabel',[],'YTickLabel',[]);
            
            hold on;
            line([0 0],get(gca,'YLim'),'Color','k');
            line(get(gca,'XLim'),[0 0],'Color','k');
            
            text(min(get(gca,'XLim')),max(get(gca,'YLim')),num2str(gridOrient(i)),'Color','r');
            
            if ismember(gridOrient(i),dat.badChans)
                set(gca,'Color','r');
            end
        end
    end
    
    %% PLOT SINGLE TRIAL RASTERS
    
    if plot_raster_flag
        figure;
        gridOrient = dat.gridOrient';
        
        for i = 1:size(gridOrient,1)*size(gridOrient,2)
            if ~isnan(gridOrient(i))
                p = plotGridPosition_new(i,size(gridOrient,1)*size(gridOrient,2),ceil(size(gridOrient,1)));
                subplot('Position',p);
                
                %         plot(taxis,squeeze(nanmean(dat.hg_stim(i,:,:),3)));
                imagesc(squeeze(dat.hg_stim(gridOrient(i),:,vertcat(trls{:})))');
                
                %         set(gca,'CLim',[-1 5],'XTickLabel',[],'YTickLabel',[]);
                set(gca,'CLim',[-1 10],'XTickLabel',[],'YTickLabel',[]);
                
                hold on;
                line([101 101],get(gca,'YLim'),'Color','k');
                
                if ismember(gridOrient(i),dat.badChans)
                    text(min(get(gca,'XLim')),max(get(gca,'YLim')),num2str(gridOrient(i)),'Color','r');
                else
                    text(min(get(gca,'XLim')),max(get(gca,'YLim')),num2str(gridOrient(i)),'Color','y');
                end
            end
        end
    end
end
