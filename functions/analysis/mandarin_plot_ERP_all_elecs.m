function mandarin_plot_ERP_all_elecs(dat,epoch,time_axis,time_lim,trls,eleclabels,array_names_unique,elec_clrs,clr)

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
