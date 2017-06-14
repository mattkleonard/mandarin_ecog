%% TO DO
%
% 3/21/17: - Pitch contour decoding
%               - Early vs. late in learning
% 4/18/17: - ERPs for correct vs. incorrect (change to fitglme?)
%          - confusion mats flipped across diagonal
%          - when do stimuli maximally diverge in terms of pitch
%          - drift diffusion model - ANOVA style
%               - which DD model parameter is affected by learning
%          - look at pitch height and direction values themselves
%          - focus learning effects on non-STG?
%          - within-block learning of tones
% 5/4/17: - Use movmean(dat.behav.accuracy(find(dat.behav.tone(plot_trls) == tones(i))),10)
%           as regressor for each tone individually (vs. binary acc)
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
elec = 182;             % which electrode you want to look at
alpha_level = 0.05/256; % statistical threshold

% stimulus parameters
blocks = [];                        % which blocks of data to use
    if strcmpi(subj,'EC137')
        blocks = [1:3 10 12:19];
    end
syls = {'bu','di','lu','ma','mi'};  % which syllables to include
spkrs = {'a','i'};                  % which speakers to include
tones = [1 2 3 4];                  % which tones to include
acc = [0 1];                        % include incorrect/correct trials

% time parameters
time_lim = [-500 1000];

% FLAGS
compare_blocks_flag = 0;            % compare 2 blocks (e.g., early vs. late)
behav_rt_analysis_flag = 0;            % run behavioral RT analysis
plot_brain_flag = 0;                % plot brain with elecs (no data)

plot_ERP_flag = 1;                  % plot ERPs on subplots
    plot_ERP_all_elecs_flag = 0;    % plot all grid electrodes
    plot_select_elecs_flag = 0;     % plot specific elec array
        elecs = {'HD'};                 % which elec array to plot (e.g., 'HD', 'AD')
    plot_single_ERP_flag = 0;           % plot one electrode ERP
    plot_single_raster_flag = 1;        % plot one electrode raster
    
regression_flag = 0;                % perform regression analysis
    run_regress_flag = 0;               % run linear regression
    feat_mat = 'height_direction_acc_int';       % which model to run ('acc','acc_trlsBlk_int',
                                        % 'accTrlsNormByBlkInt','accTrlsContInt','accStimsExpInt',
                                        % 'accTonesExpInt','height_direction','height_direction_acc_int',
                                        % 'height_direction_TonesExp_int','toneDV_RT_int','cumsumAcc')
    plot_r2_flag = 0;                   % 1=plot R2, 0=plot betas
        bs_to_plot = [6 7];                    % which betas to look at (if empty, all)
    plot_regress_all_elecs_flag = 0;    % plot regression results for all elecs
    plot_regress_brain_flag = 0;        % plot R2 or betas on brain
    plot_regress_indiv_elec_flag = 1;   % plot single elec regression results (using 'elec' above)
    depth_regression_flag = 0;          % run and plot regression for non-grid elecs

logit_regression_flag = 0;          % perform logistic regression
    run_logit_regress_flag = 0;         % run logistic regression
    plot_logit_all_elecs_flag = 1;      % plot logistic regression results for all elecs
        logit_feats_to_plot = [2];        % which features to look at
    plot_logit_brain_flag = 1;          % plot logistic regression results on brain

anova_flag = 0;                     % perform ANOVA
    run_anova_flag = 1;                 % run ANOVA
        anova_mdl = 'height_dir';          % which model to run ('syl_tone_spk','tone_block',
                                                    % 'tone_acc','syl_tone_spk_cumAcc','syl_tone_spk_expNum',
                                                    % 'height_dir_acc')
        plotConds = [];                     % which conds to plot (if [], plot all)
    plot_anova_all_elecs_flag = 1;      % plot ANOVA F-values for all elecs
    plot_anova_mean_f_flag = 1;         % plot mean F-values across elecs
    plot_anova_indiv_elec_flag = 0;     % plot single elec F-values (using 'elec' above)

classifier_flag = 0;                % perform classification
    run_classifier_flag = 1;            % run classifier (LDA)
        win_size = 10;                      % average over win_size timepoints
        kfold_xval = 10;                    % cross-validate model with kfold_xval folds
    plot_classifier_acc_flag = 1;       % plot classifier accuracy timecourses
    plot_abs_acc_flag = 0;              % plot accuracy as folds over chance

save_analysis_out_flag = 0;         % whether to save analysis outputs as .mat files

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
    if strcmpi(dat.hemi,'lh')
        x_offset = -5;
    else
        x_offset = 5;
    end
    
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
    
    % get anatomy data
    array_names = regexp(eleclabels(:,1),'[0-9]','split');
    array_names = [array_names{:,1}];
    [array_names_unique,idx] = unique(array_names(~strcmpi(array_names,'')),'stable');
    array_names = array_names(~strcmpi(array_names,''));
    elec_nums = regexp(eleclabels(:,1),['\d'],'match');
    elec_clrs = cbrewer('qual','Paired',length(array_names_unique));
    
    time_axis = dat.time_axis*1000; % [-1000:10:1000]; % [-500:10:1990];
    taxis = time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2)));

    % colors
    clr = [0.667, 0.224, 0.224 ; 0.176, 0.533, 0.176 ; 0.133, 0.40, 0.40 ; 0.667, 0.424, 0.224 ; 0.863 0.502 0.698];

end

%% GET TRIALS

[trls,plot_trls] = mandarin_get_trls(dat,trls_org,ripple_flag,compare_blocks_flag,blocks,syls,spkrs,tones,acc);

%% BEHAVIORAL ANALYSIS

if behav_rt_analysis_flag
    [r,p,all_r,all_p] = mandarin_behav_rt_analysis(dat,trls,plot_trls,clr);
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

%% PLOT ERPs

if plot_ERP_flag
    mandarin_plot_ERP(dat,epoch,time_axis,taxis,time_lim,trls,plot_trls,trls_org,elecs,blocks,syls,tones,elec,eleclabels,array_names,array_names_unique,anatomy,elec_clrs,clr,plot_ERP_all_elecs_flag,plot_select_elecs_flag,plot_single_ERP_flag,plot_single_raster_flag);
end

%% REGRESSION

if regression_flag
    
    if run_regress_flag
        b = [];
        bint = [];
        stats = [];
        sig_elecs = [];
    end

    [b,bint,stats,sig_elecs] = mandarin_regression(dat,feat_mat,taxis,time_axis,time_lim,plot_trls,syls,spkrs,tones,epoch,elec,run_regress_flag,plot_regress_all_elecs_flag,depth_regression_flag,bs_to_plot,plot_r2_flag,plot_regress_brain_flag,alpha_level,eleclabels,elec_clrs,array_names_unique,cortex,elecmatrix,x_offset,plot_regress_indiv_elec_flag,b,bint,stats,sig_elecs);
end


%% LOGISTIC REGRESSION

if logit_regression_flag
    
    if run_logit_regress_flag
        b = [];
        dev = [];
        bint = [];
        stats = [];
        coeff = [];
        t_logit = [];
        p_logit = [];
        sig_elecs = [];
    end

    [b,dev,stats,bint,coeff,t_logit,p_logit] = mandarin_logit_regression(dat,taxis,epoch,plot_trls,alpha_level,run_logit_regress_flag,plot_logit_all_elecs_flag,plot_logit_brain_flag,cortex,elecmatrix,eleclabels,x_offset,logit_feats_to_plot,b,dev,stats,bint,coeff,t_logit,p_logit,sig_elecs);
end


%% ANOVA

if anova_flag
    
    if run_anova_flag
        tbl = [];
        Fstats = [];
        pvals = [];
    end
    
    [tbl,Fstats,pvals] = mandarin_anova(dat,trls,epoch,taxis,tones,alpha_level,run_anova_flag,anova_mdl,save_analysis_out_flag,rootdir,plot_anova_all_elecs_flag,plotConds,plot_anova_mean_f_flag,plot_anova_indiv_elec_flag,tbl,Fstats,pvals);
end

%% CLASSIFIER

if classifier_flag
    [mdl,label,score,cost,mdl_acc] = mandarin_classifier(dat,trls,epoch,taxis,sigChans,win_size,kfold_xval,save_analysis_out_flag,rootdir,run_classifier_flag,plot_classifier_acc_flag,plot_abs_acc_flag,clr);
end
