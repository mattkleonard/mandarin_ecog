function [b,dev,stats,bint,coeff,t_logit,p_logit] = mandarin_logit_regression(dat,taxis,epoch,plot_trls,alpha_level,run_logit_regress_flag,plot_logit_all_elecs_flag,plot_logit_brain_flag,cortex,elecmatri,eleclabels,x_offset,logit_feats_to_plot,b,dev,stats,bint,coeff,t_logit,p_logit,sig_elecs)

bclr = [0 0 0.6 ; 0 0 0 ; 0 0.6 0 ; 0.8 0 0 ; 0.8 0 0.8; 0.6 0.6 0];

if run_logit_regress_flag
    clear b dev stats bint coeff t_logit p_logit
    
    wb = waitbar(0,'Running regression models');
    fprintf('Running regression models....\n');
    
    for i = 1:numel(dat.gridOrient)
        waitbar(i/size(dat.(epoch),1));
        for j = 1:length(taxis)
            X = [squeeze(dat.(epoch)(dat.gridOrient(i),j,plot_trls)) ...
                ];
            [b(i,j,:),dev(i,j),stats(i,j)] = glmfit(X,dat.behav.accuracy(plot_trls),...
                'binomial','link','logit');
            bint(i,j,:) = stats(i,j).se;
            coeff(i,j,:) = stats(i,j).coeffcorr(1:2);
            t_logit(i,j,:) = stats(i,j).t;
            p_logit(i,j,:) = stats(i,j).p;
        end
    end
    
end

if plot_logit_all_elecs_flag
    
    figure;
    for i = 1:numel(dat.gridOrient)
        wb = waitbar(i/size(dat.(epoch),1));
        p = plotGridPosition_new(i,numel(dat.gridOrient),ceil(sqrt(numel(dat.gridOrient))));
        ax(i) = subplot('Position',p);
        
        for k = logit_feats_to_plot
            if ~ismember(dat.gridOrient(i),dat.badChans)
                set(gca,'YLim',[min(min(min(squeeze(t_logit(:,:,logit_feats_to_plot))))) max(max(max(squeeze(t_logit(:,:,logit_feats_to_plot)))))]);
                for j = 1:size(stats,2)
                    if p_logit(dat.gridOrient(i),j,k) <= alpha_level
                        line([taxis(j) taxis(j)],get(gca,'YLim'),'Color',[0.75 0.75 0.75]);
                        hold on;
                    end
                end
                
                if ~isempty(find(p_logit(dat.gridOrient(i),:,k) <= alpha_level))
                    sig_elecs(i) = 1;
                else
                    sig_elecs(i) = 0;
                end
                
                plot(taxis,squeeze(t_logit(dat.gridOrient(i),:,k)));
                hold on;
                set(gca,'XLim',[min(taxis) max(taxis)]);
                
                elec_array = regexp(eleclabels(i,1),'[0-9]','split');
                elec_array = elec_array{:,1}(1);
                
                set(gca,'XTickLabel',[],'YTickLabel',[],...
                    'XColor',elec_clrs(find(strcmpi(elec_array,array_names_unique)),:),...
                    'YColor',elec_clrs(find(strcmpi(elec_array,array_names_unique)),:));%,...
                %             'YLim',[0 0.5]);
                
                set(gca,'YLim',[min(min(min(squeeze(t_logit(:,:,logit_feats_to_plot))))) max(max(max(squeeze(t_logit(:,:,logit_feats_to_plot)))))]);
                
            else
                set(gca,'XTickLabel',[],'YTickLabel',[],...
                    'XColor','r','YColor','r');
            end
        end
        line([0 0],get(gca,'YLim'),'Color','k');
        line(get(gca,'XLim'),[0 0],'Color','k');
        text(min(get(gca,'XLim')) + (min(get(gca,'XLim'))*0.1),...
            max(get(gca,'YLim')) - (max(get(gca,'YLim'))*0.1),...
            num2str(dat.gridOrient(i)));
    end
end
close(wb);

if plot_logit_brain_flag
    figure;
    ctmr_gauss_plot(cortex,[0 0 0],0,dat.hemi);
    hold on;
    plot_sig_elecs = find(sig_elecs);
    cmap = cbrewer('seq','Reds',ceil(max(max(squeeze(t_logit(dat.gridOrient(plot_sig_elecs),:,2)))*100)));
    for i = 1:size(t_logit,1); %length(plot_sig_elecs)
        if ~ismember(i,plot_sig_elecs)
            scatter3(elecmatrix(dat.gridOrient(i),1)+x_offset,elecmatrix(dat.gridOrient(i),2),elecmatrix(dat.gridOrient(i),3),...
                50,'k');
        else
            scatter3(elecmatrix(dat.gridOrient(i),1)+x_offset,elecmatrix(dat.gridOrient(i),2),elecmatrix(dat.gridOrient(i),3),...
                50,...
                cmap(round(max(squeeze(t_logit(dat.gridOrient(i),:,2)))*100),:),...
                'filled');
        end
    end
end
