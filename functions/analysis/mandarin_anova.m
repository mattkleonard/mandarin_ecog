function [tbl,Fstats,pvals] = mandarin_anova(dat,trls,epoch,taxis,tones,alpha_level,run_anova_flag,anova_mdl,save_analysis_out_flag,rootdir,plot_anova_all_elecs_flag,plotConds,plot_anova_mean_f_flag,plot_anova_indiv_elec_flag,tbl,Fstats,pvals)

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
    exposure_norm = (exposure - min(exposure)) / (max(exposure) - min(exposure));
    
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
                case 'tone_acc'
                    [~,tbl{i,j}] = anovan(squeeze(dat.(epoch)(i,j,plot_trls)),...
                        {dat.behav.tone(plot_trls),...
                        dat.behav.accuracy(plot_trls)},...
                        'model','full',...
                        'varnames',{'tone','acc'},...
                        'display','off');
                case 'syl_tone_spk_cumAcc'
                    cumsumTrls = [];
                    for t = 1:length(tones)
                        cumsumTrls = [cumsumTrls ; cumsum(dat.behav.accuracy(find(dat.behav.tone == tones(t))))];
                    end
                    cumsumTrls = [];
                    for t = 1:length(tones)
                        tmp_trls = find(dat.behav.tone == tones(t));
                        tmp_cumsumTrls = [];
                        for k = 1:length(tmp_trls)
                            if dat.behav.accuracy(tmp_trls(k)) == 1
                                tmp_cumsumTrls = [tmp_cumsumTrls 1];
                            else
                                tmp_cumsumTrls = [tmp_cumsumTrls -1];
                            end
                        end
                        cumsumTrls = [cumsumTrls cumsum(tmp_cumsumTrls)];
                    end
                    
                    [~,tbl{i,j}] = anovan(squeeze(dat.(epoch)(i,j,plot_trls)),...
                        {dat.behav.syllable(plot_trls),...
                        dat.behav.tone(plot_trls),...
                        dat.behav.speaker(plot_trls),...
                        cumsumTrls(plot_trls)},...
                        'model','full',...
                        'varnames',{'syllable','tone','speaker','acc'},...
                        'continuous',[4],...
                        'display','off');
                case 'syl_tone_spk_expNum'
                    [~,tbl{i,j}] = anovan(squeeze(dat.(epoch)(i,j,plot_trls)),...
                        {dat.behav.syllable(plot_trls),...
                        dat.behav.tone(plot_trls),...
                        dat.behav.speaker(plot_trls),...
                        exposure_norm(plot_trls)},...
                        'model','full',...
                        'varnames',{'syllable','tone','speaker','expNum'},...
                        'continuous',[4],...
                        'display','off');
                case 'height_dir'
                    [~,tbl{i,j}] = anovan(squeeze(dat.(epoch)(i,j,plot_trls)),...
                        {dat.behav.pitch_height(plot_trls),...
                        dat.behav.pitch_direction(plot_trls)},...
                        'model','full',...
                        'varnames',{'height','direction'},...
                        'continuous',[1 2],...
                        'display','off');
                case 'height_dir_acc'
                    [~,tbl{i,j}] = anovan(squeeze(dat.(epoch)(i,j,plot_trls)),...
                        {dat.behav.pitch_height(plot_trls),...
                        dat.behav.pitch_direction(plot_trls),...
                        dat.behav.accuracy(plot_trls)},...
                        'model','full',...
                        'varnames',{'height','direction','accuracy'},...
                        'continuous',[1 2],...
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
