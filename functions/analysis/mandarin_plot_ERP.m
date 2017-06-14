function mandarin_plot_ERP(dat,epoch,time_axis,taxis,time_lim,trls,plot_trls,trls_org,elecs,blocks,syls,tones,elec,eleclabels,array_names,array_names_unique,anatomy,elec_clrs,clr,plot_ERP_all_elecs_flag,plot_select_elecs_flag,plot_single_ERP_flag,plot_single_raster_flag)

%% PLOT ALL GRID ELECTRODES
if plot_ERP_all_elecs_flag
    fprintf('Plotting ERPs for all grid electrodes....\n');
    
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

%% PLOT SELECT ELECTRODES BY ANATOMY

if plot_select_elecs_flag
    
    fprintf('Plotting ERPs for selected electrodes....\n');
    
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
