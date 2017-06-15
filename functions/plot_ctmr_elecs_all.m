function plot_ctmr_elecs_all(rootdir,subj,hemi,elecmat,elec_num_flag,facealpha)

% subj = 'EC131';
% hemi = 'lh';
% elecmat = 'clinical';
% elec_num_flag = 1;
% facealpha = 0.1;
% 
% rootdir = '/Users/mattleonard/Documents/Research/dura/data_store2/imaging/subjects';

%%

load([rootdir '/' subj '/Meshes/' subj '_' hemi '_pial.mat']);
load([rootdir '/' subj '/elecs/' elecmat '_elecs_all.mat']);

figure('Color',[1 1 1]);
ctmr_gauss_plot(cortex,[0 0 0],0,hemi);
hold on;

array_names = regexp(eleclabels(:,1),'[0-9]','split');
array_names = [array_names{:,1}];
[array_names_unique,idx] = unique(array_names(~strcmpi(array_names,'')),'stable');
array_names = array_names(~strcmpi(array_names,''));
elec_nums = regexp(eleclabels(:,1),['\d'],'match');

if max(idx) > size(elecmatrix,1)
    idx(find(idx > size(elecmatrix,1))) = [];
end
clrs = cbrewer('qual','Paired',length(array_names_unique));
for i = 1:size(elecmatrix,1)
    elec_array = regexp(eleclabels(i,1),'[0-9]','split');
    elec_array = elec_array{:,1}(1);
    p(i) = scatter3(elecmatrix(i,1),elecmatrix(i,2),elecmatrix(i,3),50,clrs(find(strcmpi(elec_array,array_names_unique)),:),'filled');
    if elec_num_flag
        text(elecmatrix(i,1),elecmatrix(i,2),elecmatrix(i,3),char(elec_nums{i})')
    end
    hold on;
end
set(gcf,'Units','normalized');

kids = get(gca,'Children');
set(kids(end),'FaceAlpha',facealpha);

legend(p(idx),array_names(idx),'FontSize',12,'Box','off');
% leg = axes('Position',[0.85 0.6 0.1 0.4]);
