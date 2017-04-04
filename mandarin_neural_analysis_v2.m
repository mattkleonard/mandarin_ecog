%% TO DO
%
% 3/20/17: - Implement classifier over time for tone, syllable, speaker / DONE
%               - Time-locked to vowel onset/voicing onset / DONE
%          - Regression for pitch height and direction / DONE
%               - Split data into early and late within task / DONE
%          - Learning effect with exposure num instead of trial num / DONE
% 3/21/17: - Pitch contour decoding
%               - Early vs. late in learning
% 3/23/17: - Add save flag for analysis outputs / DONE
% 3/30/17: - ANOVA: Tone x Block / DONE
%               - Also do separate ANOVAs by block
%          - ANOVA: Tone x Syl x Spk x Cont_expNum / DONE

%% NOTES
% Quickly identify M/F to build relative pitch perceptual space. This
% allows tone information to not be overlapping across speakers.
%
% Is delay between syl and tone related to performance?
%
% Identify tone selective electrodes, look at behavioral and neural
% confusions over blocks.

%% PARAMETERS/FLAGS

subj = 'EC131';
local_dir_flag = 1;

if ~local_dir_flag
    rootdir = ['/Users/mattleonard/Documents/Research/pia/userdata/matt_l/mandarin/' subj];
else
    rootdir = ['/Users/mattleonard/Documents/Research/data/' subj];
end

% trial parameters
ripple_flag = 0;        % whether to run ripple task data
epoch = 'hg_stim';      % 'hg_stim','hg_vowel','hg_button','hg_feedback'
trls_org = 'tones';     % 'tones', 'syls'
elec = 170;              % which electrode you want to look at
alpha_level = 0.05/256; % statistical threshold

% stimulus parameters
blocks = [1:15];                    % which blocks of data to use
syls = {'bu','di','lu','ma','mi'};  % which syllables to include
spkrs = {'a','i'};                  % which speakers to include
tones = [1 2 3 4];                  % which tones to include
acc = [0 1];                        % whether to include correct/incorrect trials

% time parameters
time_lim = [-500 1000];

% FLAGS
behav_analysis_flag = 0;            % run behavioral RT analysis
plot_brain_flag = 0;                % plot brain with elecs (no data)
plot_ERP_flag = 0;                  % plot ERPs on subplots
compare_blocks_flag = 0;            % compare 2 blocks (e.g., early vs. late)
plot_select_elecs_flag = 0;         % plot specific elec array
    elecs = {'HD'};                     % which elec array to plot (e.g., 'HD', 'AD')
plot_single_ERP_flag = 0;           % plot one electrode ERP
plot_single_raster_flag = 0;        % plot one electrode raster
regression_flag = 0;                % perform regression analysis
    run_regress_flag = 0;               % run linear regression
    feat_mat = 7;                       % which behavioral parameters to regress
    plot_r2_flag = 0;                   % 1=plot R2, 0=plot betas
        bs_to_plot = [2:3];                 % which betas to look at
    plot_regress_all_elecs_flag = 0;    % plot regression results for all elecs
    plot_regress_brain_flag = 0;        % plot R2 or betas on brain
    plot_regress_indiv_elec_flag = 1;   % plot single elec regression results (using 'elec' above)
    depth_regression_flag = 0;          % run and plot regression for non-grid elecs
anova_flag = 1;                     % perform ANOVA
    run_anova_flag = 0;                 % run ANOVA
        anova_mdl = 'syl_tone_spk_expNum';     % which model to run ('syl_tone_spk','tone_block','syl_tone_spk_expNum')
        plotConds = [1 2 3 4 7 9 10];                     % which conds to plot (if [], plot all)
    plot_anova_all_elecs_flag = 1;     % plot ANOVA F-values for all elecs
    plot_anova_mean_f_flag = 1;         % plot mean F-values across elecs
    plot_anova_indiv_elec_flag = 0;     % plot single elec F-values (using 'elec' above)
classifier_flag = 0;                % perform classification
    run_classifier_flag = 1;            % run classifier (LDA)
        win_size = 10;                      % average over win_size timepoints
        kfold_xval = 10;                    % cross-validate model with kfold_xval folds
    plot_classifier_acc_flag = 1;       % plot classifier accuracy timecourses
    plot_abs_acc_flag = 0;              % plot accuracy as folds over chance
save_analysis_out_flag = 0;         % whether to save analysis outputs as .mat files

% colors
clr = [0.667, 0.224, 0.224 ; 0.176, 0.533, 0.176 ; 0.133, 0.40, 0.40 ; 0.667, 0.424, 0.224 ; 0.863 0.502 0.698];
% clr = [1 0.737 0.58 ; 0.435, 0.753, 0.67 ; 0.463, 0.592, 0.757 ; 1, 0.839, 0.58 ; 0.863 0.502 0.698]; % {'r','g','b','m'};
% clr = cbrewer('seq','Blues',length(blocks));

%% LOAD DATA

if ~exist('dat','var')
    fprintf('Loading data....\n');
    if ~local_dir_flag
        if ripple_flag
            load([rootdir '/data/' subj '_ripple_dat.mat']);
        else
            load([rootdir '/data/' subj '_dat.mat']);
        end
        load([rootdir '/../../../../data_store2/imaging/subjects/' subj '/elecs/TDT_elecs_all.mat']);
        load([rootdir '/../../../../data_store2/imaging/subjects/' subj '/Meshes/' subj '_' dat.hemi '_pial.mat']);
    else
        if ripple_flag
            load([rootdir '/mandarin/data/' subj '_ripple_dat.mat']);
        else
            load([rootdir '/mandarin/data/' subj '_dat.mat']);
        end
        load([rootdir '/../MRI/' subj '/elecs/TDT_elecs_all.mat']);
        load([rootdir '/../MRI/' subj '/Meshes/' subj '_' dat.hemi '_pial.mat']);
    end
end

%% GET TRIALS

trls = [];
if compare_blocks_flag
    trls{1} = find((ismember(dat.behav.block,accBlock.good{tones}))' & ...
        (ismember(dat.behav.tone,tones)) & ...
        (ismember(dat.behav.accuracy,acc)));
    trls{2} = find((ismember(dat.behav.block,accBlock.bad{tones}))' & ...
        (ismember(dat.behav.tone,tones)) & ...
        (ismember(dat.behav.accuracy,acc)));
elseif ~ripple_flag
    if strcmpi(trls_org,'tones')
        for i = 1:length(tones)
            trls{i} = find((ismember(dat.behav.block,blocks)) & ...
                (ismember(dat.behav.tone,tones(i))) & ...
                (ismember(dat.behav.accuracy,acc)));
        end
    elseif strcmpi(trls_org,'syls')
        for i = 1:length(syls)
            trls{i} = find((ismember(dat.behav.block,blocks)) & ...
                (strcmpi(dat.behav.syllable,syls{i})) & ...
                (ismember(dat.behav.accuracy,acc)));
        end
    end
else
    for i = 1:length(tones)
        trls{i} = find((ismember(dat.behav.tone,tones(i))) & ...
            (ismember(dat.behav.accuracy,acc)));
    end
end

plot_trls = [];
for i = 1:length(trls)
    plot_trls = [plot_trls ; trls{i}];
end

% get anatomy data
array_names = regexp(eleclabels(:,1),'[0-9]','split');
array_names = [array_names{:,1}];
[array_names_unique,idx] = unique(array_names(~strcmpi(array_names,'')),'stable');
array_names = array_names(~strcmpi(array_names,''));
elec_nums = regexp(eleclabels(:,1),['\d'],'match');
elec_clrs = cbrewer('qual','Paired',length(array_names_unique));

time_axis = dat.time_axis*1000; % [-1000:10:1000]; % [-500:10:1990];
taxis = time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2)));

% add pitch height and direction to dat.behav
if ~isfield(dat.behav,'pitch_height')
    if ~local_dir_flag
        [num,txt,raw] = xlsread([rootdir '/../../../../../tasks/ToneCat_ptb3/StimulusSpace.xls']);
    else
        [num,txt,raw] = xlsread([rootdir '/../../tasks/ToneCat_ptb3/StimulusSpace.xls']);
    end
    
    for i = 1:length(dat.behav.trial)
        stim = [dat.behav.syllable{i} num2str(dat.behav.tone(i)) '-' dat.behav.speaker{i} 'N.wav'];
        dat.behav.pitch_height(i) = raw{find(strcmpi(raw(:,1),['msounds/' stim])),2};
        dat.behav.pitch_direction(i) = raw{find(strcmpi(raw(:,1),['msounds/' stim])),3};
    end
    
    
    dat.behav.pitch_height = dat.behav.pitch_height';
    dat.behav.pitch_direction = dat.behav.pitch_direction';
end

% if ripple_flag
%         dat.hg_stim(find(strcmpi(array_names,'HD')),:,find(ismember(dat.behav.block,[4 5 6]))) = NaN;
%         dat.hg_button(find(strcmpi(array_names,'HD')),:,find(ismember(dat.behav.block,[4 5 6]))) = NaN;
%         dat.hg_feedback(find(strcmpi(array_names,'HD')),:,find(ismember(dat.behav.block,[4 5 6]))) = NaN;
% end

%% BEHAVIORAL ANALYSIS

if behav_analysis_flag
    figure;
    
    for i = 1:length(trls)
        scatter(dat.behav.trial(trls{i}),dat.behav.rt(trls{i}),50,clr(i,:),'filled');
        hold on;
        [r_tmp(i,:,:),p_tmp(i,:,:)] = corrcoef(dat.behav.trial(trls{i}),dat.behav.rt(trls{i}));
        r(i) = r_tmp(i,1,2);
        p(i) = p_tmp(i,1,2);
    end
    h = lsline;
    for i = 1:length(h)
        h(i).Color = clr(i,:);
        h(i).LineWidth = 2;
    end
    set(gca,'YLim',[0 7]);
    xlabel('Trial');
    ylabel('RT (sec)');
    
    [all_r,all_p] = corrcoef(dat.behav.trial(plot_trls),dat.behav.rt(plot_trls));
    fprintf('Overall RT x trial: r=%2.2g, p=%2.2g\n',all_r(2),all_p(2));
end

%% PLOT ELECS ON BRAIN

if plot_brain_flag
    fprintf('Plotting brain and elecs....\n');
    if ~local_dir_flag
        plotDir = [rootdir '/../../../../data_store2/imaging/subjects'];
    else
        plotDir = [rootdir '/../MRI'];
    end
    plot_ctmr_elecs_all(plotDir,dat.subj,dat.hemi,'TDT',1,0.2);
end

%% PLOT ALL ELECS

if plot_ERP_flag
    fprintf('Plotting ERPs....\n');
    
    clear tmp ylims;
    for i = 1:numel(dat.gridOrient)
        if ~isnan(dat.gridOrient(i))
            for j = 1:length(trls)
                tmp(i,j,:) = squeeze(nanmean(dat.(epoch)(dat.gridOrient(i),find(time_axis == time_lim(1)):find(time_axis == time_lim(2)),trls{j}),3));
            end
        end
        ylims(i,:) = [min(min(min(tmp))) max(max(max(tmp)))];
    end
    
    figure;
    wb = waitbar(0,'Plotting data');
    for i = 1:numel(dat.gridOrient)
        waitbar(i/numel(dat.gridOrient));
        p = plotGridPosition_new(i,numel(dat.gridOrient),ceil(sqrt(numel(dat.gridOrient))));
        ax(i) = subplot('Position',p);
        
        if ~ismember(dat.gridOrient(i),dat.badChans) & ~isnan(dat.gridOrient(i))
            for j = 1:length(trls)
                h(i,j) = shadedErrorBar(time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2))),...
                    squeeze(nanmean(dat.(epoch)(dat.gridOrient(i),find(time_axis == time_lim(1)):find(time_axis == time_lim(2)),trls{j}),3)),...
                    squeeze(ste(dat.(epoch)(dat.gridOrient(i),find(time_axis == time_lim(1)):find(time_axis == time_lim(2)),trls{j}),3)),...
                    {'-','Color',clr(j,:),'LineWidth',2},1);
                hold on;
            end
            axis tight;
            
            elec_array = regexp(eleclabels(i,1),'[0-9]','split');
            elec_array = elec_array{:,1}(1);
            
            set(gca,'XTickLabel',[],'YTickLabel',[],...
                'XColor',elec_clrs(find(strcmpi(elec_array,array_names_unique)),:),...
                'YColor',elec_clrs(find(strcmpi(elec_array,array_names_unique)),:),...
                'YLim',[min(min(ylims)) max(max(ylims))]);
            
            line([0 0],get(gca,'YLim'),'Color','k');
            line(get(gca,'XLim'),[0 0],'Color','k');
        else
            set(gca,'XTickLabel',[],'YTickLabel',[],...
                'XColor','r','YColor','r');
        end
        text(min(get(gca,'XLim')) + (min(get(gca,'XLim'))*0.1),...
            max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1),...
            num2str(dat.gridOrient(i)));
    end
    close(wb);
end

%% PLOT SELECT ELECTRODES

if plot_select_elecs_flag
    
    elecs_idx = [];
    for i = 1:length(elecs)
        elecs_idx{i} = find(strcmpi(elecs(i),array_names));
        for j = length(elecs_idx{i}):-1:1
            if ~cellfun(@isempty,(strfind(anatomy(elecs_idx{i}(j),4),'White-Matter')))
                elecs_idx{i}(j) = [];
            end
        end
    end
    elecs_to_plot = cell2mat(elecs_idx);
    figure;
    for i = 1:length(elecs_to_plot)
        p = plotGridPosition_new(i,length(elecs_to_plot),ceil(sqrt(length(elecs_to_plot))));
        ax(i) = subplot('Position',p);
        
        for j = 1:length(trls)
            clear badTrials
            for k = 1:length(trls{j})
                if ~isempty(find(dat.(epoch)(elecs_to_plot(i),find(time_axis == time_lim(1)):find(time_axis == time_lim(2)),trls{j}(k)) > 20))
                    badTrials(k) = 1;
                else
                    badTrials(k) = 0;
                end
            end
            
            shadedErrorBar(time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2))),...
                squeeze(nanmean(dat.(epoch)(elecs_to_plot(i),find(time_axis == time_lim(1)):find(time_axis == time_lim(2)),trls{j}(~badTrials)),3)),...
                squeeze(ste(dat.(epoch)(elecs_to_plot(i),find(time_axis == time_lim(1)):find(time_axis == time_lim(2)),trls{j}(~badTrials)),3)),...
                {'-','Color',clr(j,:),'LineWidth',2},1);
            hold on;
            %         imagesc(squeeze(dat.(epoch)(elecs(i),find(time_axis == time_lim(1)):find(time_axis == time_lim(2)),trls{j}(~badTrials)))');
        end
        line([0 0],get(gca,'YLim'),'Color','k');
        line(get(gca,'XLim'),[0 0],'Color','k');
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        
        text(min(get(gca,'XLim')) + (min(get(gca,'XLim'))*-0.1),...
            max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1),...
            anatomy(elecs_to_plot(i),1));
    end
end

%% PLOT SINGLE CHANNEL

if plot_single_ERP_flag
    figure;
    clear h;
    for j = 1:length(trls)
        h(j) = shadedErrorBar(time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2))),...
            squeeze(nanmean(dat.(epoch)(elec,find(time_axis == time_lim(1)):find(time_axis == time_lim(2)),trls{j}),3)),...
            squeeze(ste(dat.(epoch)(elec,find(time_axis == time_lim(1)):find(time_axis == time_lim(2)),trls{j}),3)),...
            {'-','Color',clr(j,:),'LineWidth',2},1);
        hold on;
    end
    axis tight;
    
    elec_array = regexp(eleclabels(elec,1),'[0-9]','split');
    elec_array = elec_array{:,1}(1);
    
    line([0 0],get(gca,'YLim'),'Color','k');
    line(get(gca,'XLim'),[0 0],'Color','k');
    
    xlabel('Time (ms)');
    ylabel('High-gamma z-score');
    title(['e' num2str(elec) ' Blocks ' num2str(min(blocks)) '-' num2str(max(blocks))]);
    
    if strcmpi(trls_org,'syls')
        legend([h.mainLine],syls);
    elseif strcmpi(trls_org,'tones')
        legend([h.mainLine],cellstr(num2str(tones(:))));
    end
end

%% PLOT SINGLE CHANNEL RASTERS

if plot_single_raster_flag
    figure;
    plot_trls = [];
    for i = 1:length(trls)
        plot_trls = [plot_trls ; trls{i}];
    end
    imagesc(squeeze(dat.(epoch)(elec,find(time_axis == time_lim(1)):find(time_axis == time_lim(2)),plot_trls))');
    colormap((cbrewer('seq','Greys',100)));
    
    set(gca,'XTick',find(taxis == 0)+(min(taxis)/10):10:find(taxis == 0)+(max(taxis)/10));
    set(gca,'XTickLabel',taxis(get(gca,'XTick')));
    hold on;
    line([find(taxis == 0) find(taxis == 0)],get(gca,'YLim'),'Color','r');
    
    for i = 1:length(plot_trls)
        if dat.behav.accuracy(plot_trls(i)) == 1
            text(max(get(gca,'XLim')),i,'*','Color',clr(dat.behav.tone(plot_trls(i)),:),'FontSize',20)
        end
    end
    
    ylabel('Trial');
    xlabel('Time (ms)');
    title(['e' num2str(elec) ' Blocks ' num2str(min(blocks)) '-' num2str(max(blocks))]);
    
%     tmp = squeeze(dat.(epoch)(elec,find(time_axis == 100):find(time_axis == 300),plot_trls))';
%     [r,p] = corrcoef(nanmean(tmp,2),dat.behav.accuracy(plot_trls));
%     fprintf('%2.2g %2.2g\n',r(2),p(2));
end

%% REGRESSION

if regression_flag
    
    bclr = [0 0 0.6 ; 0 0 0 ; 0 0.6 0 ; 0.8 0 0 ; 0.8 0 0.8];
    
    plot_trls = [];
    for i = 1:length(trls)
        plot_trls = [plot_trls ; trls{i}];
    end
        
    %%%% SET UP FEATURE MATRIX
    switch feat_mat
        case 1
            % **** main effects + learning (=acc*trial_num (by block)) ****
            X = [ones(size(plot_trls,1),1) dat.behav.accuracy(plot_trls) plot_trls dat.behav.accuracy(plot_trls).*plot_trls];
            fprintf('Feature matrix: main effects + learning (=acc*trial_num (by block))\n');
        case 2
            % **** accuracy only ****
            X = [ones(size(plot_trls,1),1) dat.behav.accuracy(plot_trls)];
            fprintf('Feature matrix: accuracy only\n');
        case 3
            % **** learning=acc*trial_num (by block, normalized) ****
            X = [ones(size(plot_trls,1),1) dat.behav.accuracy(plot_trls).*...
                (plot_trls - min(plot_trls)) / (max(plot_trls) - min(plot_trls))];
            fprintf('Feature matrix: learning=acc*trial_num (by block, normalized)\n');
        case 4
            % **** learning=acc*trial_num (continuous) ****
            X = [ones(size(plot_trls,1),1) dot(dat.behav.accuracy(plot_trls)',1:length(plot_trls),1)'];
            fprintf('Feature matrix: learning=acc*trial_num (continuous)\n');
        case 5
            % **** learning=acc*exposure_num (individual stims, continuous) ****
            exposure = zeros(1,length(dat.behav.syllable));
            unique_stims = [];
            for i = 1:length(dat.behav.syllable)
                unique_stims{i} = [dat.behav.syllable{i} '_' dat.behav.speaker{i} '_' num2str(dat.behav.tone(i))];
                exposure(i) = length(find(strcmpi(unique_stims{i},unique_stims)));
            end
            X = [ones(size(plot_trls,1),1) dot(dat.behav.accuracy(plot_trls)',exposure,1)'];
            fprintf('Feature matrix: learning=acc*exposure_num (individual stims, continuous)\n');
        case 6
            % **** learning=acc*exposure_num (tones only, continuous) ****
            exposure = zeros(1,length(dat.behav.syllable));
            unique_stims = [];
            for i = 1:length(dat.behav.syllable)
                unique_stims{i} = num2str(dat.behav.tone(i));
                exposure(i) = length(find(strcmpi(unique_stims{i},unique_stims)));
            end
            X = [ones(size(plot_trls,1),1) dot(dat.behav.accuracy(plot_trls)',exposure,1)'];
            fprintf('Feature matrix: learning=acc*exposure_num (tones only, continuous)\n');
        case 7
            % **** pitch height & pitch direction ****
            X = [ones(size(plot_trls,1),1) dat.behav.pitch_height(plot_trls) dat.behav.pitch_direction(plot_trls)];
            fprintf('Feature matrix: pitch height and pitch direction\n');
        case 8
            % **** pitch height & pitch direction & exposure num (& interactions) ****
            exposure = zeros(1,length(dat.behav.syllable));
            unique_stims = [];
            for i = 1:length(dat.behav.syllable)
                unique_stims{i} = num2str(dat.behav.tone(i));
                exposure(i) = length(find(strcmpi(unique_stims{i},unique_stims)));
            end
            %             X = [ones(size(plot_trls,1),1) dat.behav.pitch_height(plot_trls) dat.behav.pitch_direction(plot_trls) (exposure/max(exposure))'];
            X = [ones(size(plot_trls,1),1) dat.behav.pitch_height(plot_trls) dat.behav.pitch_direction(plot_trls) dot(dat.behav.pitch_height(plot_trls),(exposure/max(exposure))',2) dot(dat.behav.pitch_direction(plot_trls),(exposure/max(exposure))',2)];
            fprintf('Feature matrix: pitch height and pitch direction & exposure num (& interactions)\n');
        case 9
            % **** tone & RT (& interaction) ****
            X = [ones(size(plot_trls,1),1) dummyvar(dat.behav.tone(plot_trls)).*dat.behav.rt(plot_trls)];
            fprintf('Feature matrix: Tone (dummyvar) and RT interaction (no main effects)\n');

    end
    %%%% END SET UP FEATURE MATRIX
        
    if run_regress_flag
        clear b bint stats
        
        wb = waitbar(0,'Running regression models');
        fprintf('Running regression models....\n');

        for i = 1:numel(dat.gridOrient)
            waitbar(i/size(dat.(epoch),1));
            for j = 1:length(taxis)
                [b(i,j,:),bint(i,j,:,:),~,~,stats(i,j,:)] = regress(squeeze(dat.(epoch)(dat.gridOrient(i),j,plot_trls)),X);
            end
        end
    end
    
    if plot_regress_all_elecs_flag
        figure;
        for i = 1:numel(dat.gridOrient)
            waitbar(i/size(dat.(epoch),1));
            p = plotGridPosition_new(i,numel(dat.gridOrient),ceil(sqrt(numel(dat.gridOrient))));
            ax(i) = subplot('Position',p);
            
            if ~ismember(dat.gridOrient(i),dat.badChans)
                if plot_r2_flag
                    set(gca,'YLim',[0 0.35]);
                else
                    set(gca,'YLim',[min(min(min(squeeze(b(:,:,bs_to_plot))))) max(max(max(squeeze(b(:,:,bs_to_plot)))))]);
                    %                 set(gca,'YLim',[-0.01 0.01]);
                end
                for j = 1:size(stats,2)
                    if stats(dat.gridOrient(i),j,3) <= alpha_level
                        line([taxis(j) taxis(j)],get(gca,'YLim'),'Color',[0.75 0.75 0.75]);
                        hold on;
                    end
                end
                
                if ~isempty(find(stats(dat.gridOrient(i),:,3) <= alpha_level))
                    sig_elecs(i) = 1;
                else
                    sig_elecs(i) = 0;
                end
                
                if plot_r2_flag
                    plot(taxis,squeeze(stats(dat.gridOrient(i),:,1)));
                    hold on;
                else
                    for j = bs_to_plot % 4:4; %size(b,3)
                        %                     shadedErrorBar(time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2))),...
                        %                         squeeze(b(i,:,j)),...
                        %                         [squeeze(bint(i,:,j,1)) ; squeeze(bint(i,:,j,2))],...
                        %                         {'-','Color',clr(j,:),'LineWidth',2},1);
                        [hl(i,j),hp(i,j)] = errorarea_glm(time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2))),...
                            squeeze(b(dat.gridOrient(i),:,j)),...
                            squeeze(bint(dat.gridOrient(i),:,j,1)),squeeze(bint(dat.gridOrient(i),:,j,2)));
                        set(hl(i,j),'Color',bclr(j,:),'LineWidth',2);
                        set(hp(i,j),'FaceColor',bclr(j,:));
                        hold on;
                    end
                end
                %         axis tight;
                set(gca,'XLim',[min(taxis) max(taxis)]);
                
                elec_array = regexp(eleclabels(i,1),'[0-9]','split');
                elec_array = elec_array{:,1}(1);
                
                set(gca,'XTickLabel',[],'YTickLabel',[],...
                    'XColor',elec_clrs(find(strcmpi(elec_array,array_names_unique)),:),...
                    'YColor',elec_clrs(find(strcmpi(elec_array,array_names_unique)),:));%,...
                %             'YLim',[0 0.5]);
                
                if plot_r2_flag
                    set(gca,'YLim',[0 0.35]);
                else
                    set(gca,'YLim',[min(min(min(squeeze(b(:,:,bs_to_plot))))) max(max(max(squeeze(b(:,:,bs_to_plot)))))]);
                    %                 set(gca,'YLim',[-0.01 0.01]);
                end
                line([0 0],get(gca,'YLim'),'Color','k');
                line(get(gca,'XLim'),[0 0],'Color','k');
            else
                set(gca,'XTickLabel',[],'YTickLabel',[],...
                    'XColor','r','YColor','r');
            end
            text(min(get(gca,'XLim')) + (min(get(gca,'XLim'))*0.1),...
                max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1),...
                num2str(dat.gridOrient(i)));
        end
        close(wb);
    end
    
    fprintf('Max R^2 value: %2.2g\n',max(max(squeeze(stats(:,:,1)))));
    
    %%
    if plot_regress_brain_flag
        figure;
        ctmr_gauss_plot(cortex,[0 0 0],0,dat.hemi);
        hold on;
        plot_sig_elecs = find(sig_elecs);
        if plot_r2_flag
            cmap = cbrewer('seq','Reds',ceil(max(max(squeeze(stats(dat.gridOrient(plot_sig_elecs),:,1)))*100)));
            for i = 1:size(stats,1); %length(plot_sig_elecs)
                if ~ismember(i,plot_sig_elecs)
                    scatter3(elecmatrix(dat.gridOrient(i),1)-5,elecmatrix(dat.gridOrient(i),2),elecmatrix(dat.gridOrient(i),3),...
                        50,'k');
                else
                    plot_dat_elec = max(squeeze(stats(dat.gridOrient(i),:,1)));
                    scatter3(elecmatrix(dat.gridOrient(i),1)-5,elecmatrix(dat.gridOrient(i),2),elecmatrix(dat.gridOrient(i),3),...
                        50,...
                        cmap(round(max(squeeze(stats(dat.gridOrient(i),:,1)))*100),:),...
                        'filled');
                end
            end
        elseif ~plot_r2_flag && length(bs_to_plot) == 1
            for i = 1:size(b,1)
                if ismember(i,plot_sig_elecs)
                    [~,tmp_idx] = max(abs(squeeze(b(dat.gridOrient(i),:,bs_to_plot))));
                    plot_dat(i) = squeeze(b(dat.gridOrient(i),tmp_idx,bs_to_plot));
                else
                    plot_dat(i) = NaN;
                end
            end
            plot_dat = (plot_dat - min(plot_dat)) / (max(plot_dat) - min(plot_dat));
            cmap = flipud(cbrewer('div','RdBu',101)); % ceil(max(max(squeeze(b(dat.gridOrient(plot_sig_elecs),:,bs_to_plot)))*100))));
            for i = 1:size(b,1); %length(plot_sig_elecs)
                if ~ismember(i,plot_sig_elecs)
                    scatter3(elecmatrix(dat.gridOrient(i),1)-5,elecmatrix(dat.gridOrient(i),2),elecmatrix(dat.gridOrient(i),3),...
                        50,'k');
                else
                    plot_dat_elec = plot_dat(i); % squeeze(b(dat.gridOrient(i),tmp_idx,bs_to_plot));
                    scatter3(elecmatrix(dat.gridOrient(i),1)-5,elecmatrix(dat.gridOrient(i),2),elecmatrix(dat.gridOrient(i),3),...
                        50,...
                        cmap(round(plot_dat_elec*100)+1,:),...
                        'filled');
                end
            end
            
        end
        title(epoch);
    end
    
    if plot_regress_indiv_elec_flag
        figure;
        
        for j = bs_to_plot % 4:4; %size(b,3)
            [hl(i,j),hp(i,j)] = errorarea_glm(time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2))),...
                squeeze(b(elec,:,j)),...
                squeeze(bint(elec,:,j,1)),squeeze(bint(elec,:,j,2)));
            set(hl(i,j),'Color',bclr(j,:),'LineWidth',2);
            set(hp(i,j),'FaceColor',bclr(j,:));
            hold on;
            
            axis tight;
            line([0 0],get(gca,'YLim'),'Color','k');
            line(get(gca,'XLim'),[0 0],'Color','k');
            set(gca,'XLim',[min(taxis) max(taxis)]);
            
            xlabel('Time (ms)');
            ylabel('Beta');
            
            title(['ch' num2str(elec)]);
        end
    end
end

%% FOR DEPTH ELECTRODE REGRESSION

if depth_regression_flag
    
    bs_to_plot = [2];
    bclr = [0 0 0.6 ; 0 0 0 ; 0 0.6 0 ; 0.8 0 0];
    
    plot_trls = [];
    for i = 1:length(trls)
        plot_trls = [plot_trls ; trls{i}];
    end
    
    %     X = [ones(size(plot_trls,1),1) dat.behav.accuracy(plot_trls) plot_trls dat.behav.accuracy(plot_trls).*plot_trls];
    %     X = [ones(size(plot_trls,1),1) dat.behav.accuracy(plot_trls)];
    X = [ones(size(plot_trls,1),1) dat.behav.accuracy(plot_trls).*...
        (plot_trls - min(plot_trls)) / (max(plot_trls) - min(plot_trls))];
    figure;
    wb = waitbar(0,'Plotting data');
    fprintf('Running regression models....\n');
    
    if run_regress_flag
        clear b bint stats
        
        for i = 257:320; % 1:numel(dat.gridOrient)
            for j = 1:length(taxis)
                [b(i-256,j,:),bint(i-256,j,:,:),~,~,stats(i-256,j,:)] = regress(squeeze(dat.(epoch)(i,j,plot_trls)),X);
            end
        end
    end
    for i = 1:size(stats,1); % 1:numel(dat.gridOrient)
        waitbar(i/size(dat.(epoch),1));
        p = plotGridPosition_new(i,size(stats,1),ceil(sqrt(size(stats,1))));
        ax(i) = subplot('Position',p);
        
        %         if ~ismember(i,dat.badChans)
        if plot_r2_flag
            set(gca,'YLim',[0 0.35]);
        else
            set(gca,'YLim',[min(min(min(squeeze(b(:,:,bs_to_plot))))) max(max(max(squeeze(b(:,:,bs_to_plot)))))]);
            %                 set(gca,'YLim',[-0.01 0.01]);
        end
        for j = 1:size(stats,2)
            if stats(i,j,3) <= alpha_level
                line([taxis(j) taxis(j)],get(gca,'YLim'),'Color',[0.75 0.75 0.75]);
                hold on;
            end
        end
        
        if ~isempty(find(stats(i,:,3) <= alpha_level))
            sig_elecs(i) = 1;
        else
            sig_elecs(i) = 0;
        end
        
        if plot_r2_flag
            plot(taxis,squeeze(stats(i,:,1)));
            hold on;
        else
            for j = bs_to_plot % 4:4; %size(b,3)
                %                     shadedErrorBar(time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2))),...
                %                         squeeze(b(i,:,j)),...
                %                         [squeeze(bint(i,:,j,1)) ; squeeze(bint(i,:,j,2))],...
                %                         {'-','Color',clr(j,:),'LineWidth',2},1);
                [hl(i,j),hp(i,j)] = errorarea_glm(time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2))),...
                    squeeze(b(i,:,j)),...
                    squeeze(bint(i,:,j,1)),squeeze(bint(i,:,j,2)));
                set(hl(i,j),'Color',bclr(j,:),'LineWidth',2);
                set(hp(i,j),'FaceColor',bclr(j,:));
                hold on;
            end
        end
        %         axis tight;
        set(gca,'XLim',[min(taxis) max(taxis)]);
        
        elec_array = regexp(eleclabels(i+256,1),'[0-9]','split');
        elec_array = elec_array{:,1}(1);
        
        set(gca,'XTickLabel',[],'YTickLabel',[],...
            'XColor',elec_clrs(find(strcmpi(elec_array,array_names_unique)),:),...
            'YColor',elec_clrs(find(strcmpi(elec_array,array_names_unique)),:));%,...
        %             'YLim',[0 0.5]);
        
        if plot_r2_flag
            set(gca,'YLim',[0 0.35]);
        else
            set(gca,'YLim',[min(min(min(squeeze(b(:,:,bs_to_plot))))) max(max(max(squeeze(b(:,:,bs_to_plot)))))]);
            %                 set(gca,'YLim',[-0.01 0.01]);
        end
        line([0 0],get(gca,'YLim'),'Color','k');
        line(get(gca,'XLim'),[0 0],'Color','k');
        %         else
        %             set(gca,'XTickLabel',[],'YTickLabel',[],...
        %                 'XColor','r','YColor','r');
        %         end
        text(min(get(gca,'XLim')) + (min(get(gca,'XLim'))*-0.1),...
            max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1),...
            anatomy(i+256,1));
    end
    close(wb);
end


%% ANOVA

if anova_flag
    
    plot_trls = [];
    for i = 1:length(trls)
        plot_trls = [plot_trls ; trls{i}];
    end
    
    % run ANOVA
    if run_anova_flag
        fprintf('Running ANOVA....\n');
        
        clear p tbl Fstats pvals;
        tbl = cell(size(dat.(epoch),1),length(taxis));
        exposure = zeros(1,length(dat.behav.syllable));
        unique_stims = [];
        for i = 1:length(dat.behav.syllable)
            unique_stims{i} = num2str(dat.behav.tone(i));
            exposure(i) = length(find(strcmpi(unique_stims{i},unique_stims)));
        end
        
        for i = 1:size(dat.(epoch),1)
            fprintf('Channel [%d] of [%d]\n',i,size(dat.(epoch),1));
            textprogressbar(' ');
            for j = 1:length(taxis)
                textprogressbar((j/length(taxis))*100);
                switch anova_mdl
                    case 'syl_tone_spk'
                        [~,tbl{i,j}] = anovan(squeeze(dat.(epoch)(i,j,plot_trls)),...
                            {dat.behav.syllable(plot_trls),...
                            dat.behav.tone(plot_trls),...
                            dat.behav.speaker(plot_trls)},...
                            'model','full',...
                            'varnames',{'syllable','tone','speaker'},...
                            'display','off');
                    case 'tone_block'
                        [~,tbl{i,j}] = anovan(squeeze(dat.(epoch)(i,j,plot_trls)),...
                            {dat.behav.tone(plot_trls),...
                            dat.behav.block(plot_trls)},...
                            'model','full',...
                            'varnames',{'tone','block'},...
                            'display','off');
                    case 'syl_tone_spk_expNum'
                        [~,tbl{i,j}] = anovan(squeeze(dat.(epoch)(i,j,plot_trls)),...
                            {dat.behav.syllable(plot_trls),...
                            dat.behav.tone(plot_trls),...
                            dat.behav.speaker(plot_trls),...
                            exposure(plot_trls)},...
                            'model','full',...
                            'varnames',{'syllable','tone','speaker','expNum'},...
                            'continuous',[4],...
                            'display','off');
                end
                Fstats(i,j,:) = cell2mat(tbl{i,j}(2:end-2,6));
                pvals(i,j,:) = cell2mat(tbl{i,j}(2:end-2,7));
            end
            textprogressbar('Done.');
        end
        % get sigChans from ANOVA
        sigChans = [];
        for i = 1:size(pvals,1)
            if ~ismember(i,dat.badChans) & any(find(squeeze(pvals(i,:,:)) <= alpha_level))
                sigChans = [sigChans i];
            end
        end

        if save_analysis_out_flag & ~exist([rootdir '/mandarin/data/anova/' subj '_ANOVA_' anova_mdl '.mat'])
            save([rootdir '/mandarin/data/anova/' subj '_ANOVA_' anova_mdl '.mat'],'tbl','Fstats','pvals','-v7.3');
            save([rootdir '/mandarin/data/anova/' subj '_sigChans_' anova_mdl '.mat'],'sigChans','-v7.3');
        end

    end
    
    % plot F-values
    if plot_anova_all_elecs_flag
        
        varnames = tbl{1}(:,:,1);
        if isempty(plotConds)
            plotConds = 2:size(varnames,1)-2;
            varnames = varnames(plotConds,1);
            plotConds = plotConds - 1;
        else
            varnames = varnames(plotConds+1,1);
        end
                
        cmap = cbrewer('qual','Set1',size(Fstats,3));
        
        figure;
        for i = 1:numel(dat.gridOrient)
            p1 = plotGridPosition_new(i,numel(dat.gridOrient),ceil(sqrt(numel(dat.gridOrient))));
            subplot('Position',p1);
            
            if ~ismember(dat.gridOrient(i),dat.badChans)
                for j = plotConds
                    plot(taxis,Fstats(dat.gridOrient(i),:,j),'Color',cmap(find(plotConds == j),:));
                    hold on;
                    
                    for k = 1:length(taxis)
                        if pvals(dat.gridOrient(i),k,j) <= alpha_level
                            scatter(taxis(k),Fstats(dat.gridOrient(i),k,j),20,cmap(find(plotConds == j),:),'filled');
%                             scatter(taxis(k),Fstats(dat.gridOrient(i),k,j),20,'k','filled');
                        end
                    end
                end
                
                axis tight;
                
                set(gca,'YLim',[min(min(min(Fstats(~ismember(1:size(dat.(epoch),1),dat.badChans),:,plotConds)))) max(max(max(Fstats(~ismember(1:size(dat.(epoch),1),dat.badChans),:,plotConds))))]);
                set(gca,'XTickLabel',[],'YTickLabel',[]);
                
                line([0 0],get(gca,'YLim'),'Color','k');
                
                text(0,max(get(gca,'YLim')) - (max(get(gca,'YLim')) * 0.1),num2str(dat.gridOrient(i)));
            else
                axis off;
            end
        end
    end
    
    
    % mean F-values
    if plot_anova_mean_f_flag
        figure;
        for i = plotConds % 1:size(Fstats,3)
            h(i) = shadedErrorBar(taxis,squeeze(mean(Fstats(:,:,i),1)),...
                squeeze(ste(Fstats(:,:,i),1)),...
                {'Color',cmap(find(plotConds == i),:),'LineWidth',2},1);
            hold on;
        end
        
        axis tight;
        line([0 0],get(gca,'YLim'),'Color','k');
        title('Mean F-values');
        ylabel('F-value');
        xlabel('Time (samples)');
        
        legend([h(plotConds).mainLine],varnames);
    end
    
    % plot individual electrodes
    if plot_anova_indiv_elec_flag
        
        varnames = tbl{1}(:,:,1);
        if isempty(plotConds)
            plotConds = 2:size(varnames,1)-2;
            varnames = varnames(plotConds,1);
            plotConds = plotConds - 1;
        else
            varnames = varnames(plotConds,1);
        end

        clear h;
        
        figure;
        for i = 1:size(Fstats,3)
            h(i) = plot(taxis,squeeze(Fstats(elec,:,i)),...
                'Color',cmap(plotConds(i),:),'LineWidth',2);
            hold on;
            for k = 1:length(taxis)
                if pvals(elec,k,i) <= alpha_level
                    scatter(taxis(k),Fstats(elec,k,i),50,cmap(i,:),'filled');
                end
            end
        end
        
        axis tight;
        line([0 0],get(gca,'YLim'),'Color','k');
        title(['ch' num2str(elec)]);
        ylabel('F-value');
        xlabel('Time (samples)');
        
        legend([h],varnames)
    end
end

%% CLASSIFIER

if classifier_flag
    
    if run_classifier_flag
        
%         if ~exist([rootdir '/mandarin/data/anova/' subj '_sigChans.mat'])
%             error('Please run ANOVA to find sigChans first');
%         else
%             load([rootdir '/mandarin/data/anova/' subj '_sigChans.mat']);
%         end
                
        clear X Y mdl label score cost mdl_acc
        
        plot_trls = [];
        for i = 1:length(trls)
            plot_trls = [plot_trls ; trls{i}];
        end

        % remove NaNs
        inputDat = dat.(epoch)(:,:,plot_trls);
        Y.syllable = dat.behav.syllable(plot_trls);
        Y.speaker = dat.behav.speaker(plot_trls);
        Y.tone = regexp(sprintf('%i ',dat.behav.tone(plot_trls)),'(\d+)','match')';
        flds = fieldnames(Y);
        for i = size(inputDat,3):-1:1
            if any(any(isnan(dat.(epoch)(sigChans,:,i))))
                inputDat(:,:,i) = [];
                for j = 1:length(flds)
                    Y.(flds{j})(i,:) = [];
                end
            end
        end
        
        % run classifier
        for i = win_size+1:size(dat.(epoch),2)-win_size-1
            fprintf('Computing LDA model for timepoint [%d] of [%d]\n',i-win_size,size(dat.(epoch),2)-(win_size*2));
            
            X = squeeze(nanmean(inputDat(sigChans,i-win_size:i,:),2))';
            
            for j = 1:length(flds)
                mdl{j,i-win_size} = fitcdiscr(X,Y.(flds{j}),'Prior','empirical','KFold',kfold_xval);
                
                [label(:,i-win_size,j),score{j}(:,:,i-win_size),cost{j}(:,:,i-win_size)] = kfoldPredict(mdl{j,i-win_size});
                
                mdl_acc(i-win_size,j) = length(find(cellfun(@strcmpi,Y.(flds{j}),label(:,i-win_size,j)))) / size(Y.(flds{j}),1);
            end
        end
        
        if save_analysis_out_flag & ~exist([rootdir '/mandarin/data/classifier/' subj '_classifier_syllable_speaker_tone_' epoch '.mat'])
            save([rootdir '/mandarin/data/classifier/' subj '_classifier_syllable_speaker_tone_' epoch '.mat'],'mdl','label','score','cost','mdl_acc','-v7.3');
        end
    end
    
    %% PLOT CLASSIFIER OUTPUT
    
    if plot_classifier_acc_flag
        figure;
        
        if plot_abs_acc_flag
            for i = 1:size(mdl_acc,2)
                plot(taxis(win_size+1:size(dat.(epoch),2)-win_size-1),mdl_acc(1:length(taxis(win_size:end)),i),'Color',clr(i,:),'LineWidth',2);
                hold on;
            end
            axis tight;
            for i = 1:size(mdl_acc,2)
                line(get(gca,'XLim'),[1/length(unique(Y.(flds{i}))) 1/length(unique(Y.(flds{i})))],'Color',clr(i,:),'LineStyle','--');
            end
            
            ylabel('Classification Accuracy');
        else
            for i = 1:size(mdl_acc,2)
                plot(taxis(win_size+1:size(dat.(epoch),2)-win_size-1),(mdl_acc(:,i) / (1/length(unique(Y.(flds{i}))))),...
                    'Color',clr(i,:),'LineWidth',2);
                hold on;
            end
            axis tight;
            line(get(gca,'XLim'),[1 1],'Color','k');
            
            ylabel('Accuracy over chance');
        end
        
        line([0 0],get(gca,'YLim'),'Color','k');
        
        xlabel('Time (ms)');
        
        legend(flds);
    end
end
