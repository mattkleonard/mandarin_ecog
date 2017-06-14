function [mdl,label,score,cost,mdl_acc] = mandarin_classifier(dat,trls,epoch,taxis,sigChans,win_size,kfold_xval,save_analysis_out_flag,rootdir,run_classifier_flag,plot_classifier_acc_flag,plot_abs_acc_flag,clr)

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
