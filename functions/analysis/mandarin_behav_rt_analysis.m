function [r,p,all_r,all_p] = mandarin_behav_rt_analysis(dat,trls,plot_trls,clr)

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
