function [b,bint,stats,sig_elecs] = mandarin_regression(dat,feat_mat,taxis,time_axis,time_lim,plot_trls,syls,spkrs,tones,epoch,elec,run_regress_flag,plot_regress_all_elecs_flag,depth_regression_flag,bs_to_plot,plot_r2_flag,plot_regress_brain_flag,alpha_level,eleclabels,elec_clrs,array_names_unique,cortex,elecmatrix,x_offset,plot_regress_indiv_elec_flag,b,bint,stats,sig_elecs)

bclr = [0 0 0.6 ; 0 0 0 ; 0 0.6 0 ; 0.8 0 0 ; 0.8 0 0.8; 0.6 0.6 0 ; 0 0.6 0.6];

%%%% SET UP FEATURE MATRIX
switch feat_mat
    case 'acc'
        X = [ones(size(plot_trls,1),1) dat.behav.accuracy(plot_trls)];
        feats = {'acc'};
        fprintf('Feature matrix: acc\n');
    case 'acc_trlsBlk_int'
        X = [ones(size(plot_trls,1),1) dat.behav.accuracy(plot_trls) plot_trls dat.behav.accuracy(plot_trls).*plot_trls];
        feats = {'acc','trlsBlk','acc*trlsBlk'};
        fprintf('Feature matrix: acc_trlsBlk_int\n');
    case 'accTrlsNormByBlkInt'
        X = [ones(size(plot_trls,1),1) dat.behav.accuracy(plot_trls).*...
            (plot_trls - min(plot_trls)) / (max(plot_trls) - min(plot_trls))];
        feats = {'acc*trlsNormByBlk'};
        fprintf('Feature matrix: accTrlsNormByBlkInt\n');
    case 'accTrlsContInt'
        X = [ones(size(plot_trls,1),1) dot(dat.behav.accuracy(plot_trls)',1:length(plot_trls),1)'];
        feats = {'acc*trlsCont'};
        fprintf('Feature matrix: accTrlsContInt\n');
    case 'accStimsExpInt'
        exposure = zeros(1,length(dat.behav.syllable(plot_trls)));
        unique_stims = [];
        for i = 1:length(dat.behav.syllable(plot_trls))
            unique_stims{i} = num2str(dat.behav.tone(plot_trls(i)));
            exposure(i) = length(find(strcmpi(unique_stims{i},unique_stims)));
        end
        X = [ones(size(plot_trls,1),1) dot(dat.behav.accuracy(plot_trls)',exposure,1)'];
        feats = {'acc*stimsExp'};
        fprintf('Feature matrix: accStimsExpInt\n');
    case 'accTonesExpInt'
        exposure = zeros(1,length(dat.behav.syllable(plot_trls)));
        unique_stims = [];
        for i = 1:length(dat.behav.syllable(plot_trls))
            unique_stims{i} = num2str(dat.behav.tone(plot_trls(i)));
            exposure(i) = length(find(strcmpi(unique_stims{i},unique_stims)));
        end
        X = [ones(size(plot_trls,1),1) dot(dat.behav.accuracy(plot_trls)',exposure,1)'];
        feats = {'acc*tonesExp'};
        fprintf('Feature matrix: accTonesExpInt\n');
    case 'height_direction'
        X = [ones(size(plot_trls,1),1) dat.behav.pitch_height(plot_trls) dat.behav.pitch_direction(plot_trls) dot(dat.behav.pitch_height(plot_trls),dat.behav.pitch_direction(plot_trls),2)];
        feats = {'height','direction','height*direction'};
        fprintf('Feature matrix: height_direction\n');
    case 'height_direction_acc_int'
        X = [ones(size(plot_trls,1),1) dat.behav.pitch_height(plot_trls) dat.behav.pitch_direction(plot_trls) dat.behav.accuracy(plot_trls) dot(dat.behav.pitch_height(plot_trls),dat.behav.pitch_direction(plot_trls),2) dot(dat.behav.pitch_height(plot_trls),dat.behav.accuracy(plot_trls),2) dot(dat.behav.pitch_direction(plot_trls),dat.behav.accuracy(plot_trls),2)];
        feats = {'height','direction','acc','height*direction','height*acc','direction*acc'};
        %         X = [ones(size(plot_trls,1),1) dot(dat.behav.pitch_height(plot_trls),dat.behav.accuracy(plot_trls),2) dot(dat.behav.pitch_direction(plot_trls),dat.behav.accuracy(plot_trls),2)];
        feats = {'height*acc','direction*acc'};
        fprintf('Feature matrix: height_direction_acc\n');
    case 'height_direction_TonesExp_int'
        exposure = zeros(1,length(dat.behav.syllable(plot_trls)));
        unique_stims = [];
        for i = 1:length(dat.behav.syllable(plot_trls))
            unique_stims{i} = num2str(dat.behav.tone(plot_trls(i)));
            exposure(i) = length(find(strcmpi(unique_stims{i},unique_stims)));
        end
        %             X = [ones(size(plot_trls,1),1) dat.behav.pitch_height(plot_trls) dat.behav.pitch_direction(plot_trls) (exposure/max(exposure))'];
        X = [ones(size(plot_trls,1),1) dat.behav.pitch_height(plot_trls) dat.behav.pitch_direction(plot_trls) dot(dat.behav.pitch_height(plot_trls),(exposure/max(exposure))',2) dot(dat.behav.pitch_direction(plot_trls),(exposure/max(exposure))',2)];
        feats = {'height','direction','height*tonesExp','direction*tonesExp'};
        fprintf('Feature matrix: height_direction_TonesExp_int\n');
    case 'toneDV_RT_int'
        X = [ones(size(plot_trls,1),1) dummyvar(dat.behav.tone(plot_trls)).*dat.behav.rt(plot_trls)];
        feats = {'toneDV','RT','toneDV*RT'};
        fprintf('Feature matrix: toneDV_RT_int\n');
    case 'cumsumAcc'
        cumsumTrls = [];
        for t = 1:length(tones)
            cumsumTrls = [cumsumTrls ; cumsum(dat.behav.accuracy(find(dat.behav.tone == tones(t))))];
        end
        X = [ones(size(plot_trls,1),1) dummyvar(dat.behav.tone(plot_trls)).*cumsumTrls(plot_trls)];
        feats = {'toneDV*cumsumAcc'};
        fprintf('Feature matrix: cumsumAcc\n');
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
    close(wb);
end

if isempty(sig_elecs)
    for i = 1:numel(dat.gridOrient)
        if ~isempty(find(stats(dat.gridOrient(i),:,3) <= alpha_level))
            sig_elecs(i) = 1;
        else
            sig_elecs(i) = 0;
        end
    end
end
fprintf('Max R^2 value: %2.2g\n',max(max(squeeze(stats(:,:,1)))));

if isempty(bs_to_plot)
    bs_to_plot = 2:size(b,3);
end

if plot_regress_all_elecs_flag
    figure;
    for i = 1:numel(dat.gridOrient)
        wb = waitbar(i/size(dat.(epoch),1));
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
                set(gca,'YLim',[0 max(max(squeeze(stats(:,:,1))))]);
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

%% PLOT R2 OR WEIGHTS ON BRAIN

if plot_regress_brain_flag
    figure;
    ctmr_gauss_plot(cortex,[0 0 0],0,dat.hemi);
    hold on;
    plot_sig_elecs = find(sig_elecs);
    if plot_r2_flag
        cmap = cbrewer('seq','Reds',ceil(max(max(squeeze(stats(dat.gridOrient(plot_sig_elecs),:,1)))*100)));
        for i = 1:size(stats,1); %length(plot_sig_elecs)
            if ~ismember(i,plot_sig_elecs)
                scatter3(elecmatrix(dat.gridOrient(i),1)+x_offset,elecmatrix(dat.gridOrient(i),2),elecmatrix(dat.gridOrient(i),3),...
                    50,'k');
            else
                plot_dat_elec = max(squeeze(stats(dat.gridOrient(i),:,1)));
                scatter3(elecmatrix(dat.gridOrient(i),1)+x_offset,elecmatrix(dat.gridOrient(i),2),elecmatrix(dat.gridOrient(i),3),...
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
                scatter3(elecmatrix(dat.gridOrient(i),1)+x_offset,elecmatrix(dat.gridOrient(i),2),elecmatrix(dat.gridOrient(i),3),...
                    50,'k');
            else
                plot_dat_elec = plot_dat(i); % squeeze(b(dat.gridOrient(i),tmp_idx,bs_to_plot));
                scatter3(elecmatrix(dat.gridOrient(i),1)+x_offset,elecmatrix(dat.gridOrient(i),2),elecmatrix(dat.gridOrient(i),3),...
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
        [hl(j),hp(j)] = errorarea_glm(time_axis(find(time_axis == time_lim(1)):find(time_axis == time_lim(2))),...
            squeeze(b(elec,:,j)),...
            squeeze(bint(elec,:,j,1)),squeeze(bint(elec,:,j,2)));
        set(hl(j),'Color',bclr(j,:),'LineWidth',2);
        set(hp(j),'FaceColor',bclr(j,:));
        hold on;
        
        axis tight;
        line([0 0],get(gca,'YLim'),'Color','k');
        line(get(gca,'XLim'),[0 0],'Color','k');
        set(gca,'XLim',[min(taxis) max(taxis)]);
        
        xlabel('Time (ms)');
        ylabel('Beta');
        
        title(['ch' num2str(elec)]);
    end
    legend([hl(bs_to_plot)],feats(bs_to_plot-1));
end

%% FOR DEPTH ELECTRODE REGRESSION

if depth_regression_flag
    
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

end
