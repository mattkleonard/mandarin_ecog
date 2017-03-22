subj = 'CH';
hemi = 'lh';
plot_elecs_flag = 0;
plot_dat_flag = 0;
roi = {'superiorfrontal','caudalmiddlefrontal','parstriangularis','parsopercularis','superiortemporal','precentral','postcentral','superiorparietal','supramarginal'};
clr = [0 1 0.9];
% clrs = [-0.5 -0.5 -0.5];
clrs = linspace(-1,1,length(roi));
elecsize = 50;
rootdir = '/Users/mattleonard/Documents/Research/data/MRI';

load([rootdir '/' subj '/Meshes/' subj '_' hemi '_pial.mat']);

ctab_fid = fopen([rootdir '/aparc.annot.ctab']);
ctab = textscan(ctab_fid,'%d%s%d%d%d%d');
roi_names = ctab{2};
% roi = roi_names(2:end);
% clrs = linspace(-1,1,length(roi));

rndColorOrder = randperm(length(clrs));
clrs = clrs(rndColorOrder);
clrs(find(clrs == 0)) = 0.25;

fid = fopen([rootdir '/' subj '/' hemi '.aparc.annot.dpv']);
vert = textscan(fid,'%d%d%d%d%d');
vert = vert{5};

figure;
ctmr_gauss_plot(cortex,[0 0 0],0,hemi);

% cmap = [[0.6 0.6 0.6] ; cbrewer('qual','Set1',length(roi))];
% colormap(flipud(cmap));

kids = get(gca,'Children');

for i = 1:length(roi)
    roi_idx = find(strcmpi(roi_names,roi{i}));
    
    roi_verts{i} = find(vert == roi_idx);
    kids(2).FaceVertexCData(roi_verts{i},:) = repmat(clrs(i),length(roi_verts{i}),1);
%     kids(2).FaceVertexCData(roi_verts{i},1:3) = repmat([0.2 0.8 0.5],length(roi_verts{i}),1);
end

% kids(2).FaceAlpha = 'flat';
% kids(2).FaceVertexAlphaData(roi_verts) = 0.5;
% colormap([[0.9412    0.9412    0.9412] ; clr]);
% caxis([0 1]);

if plot_elecs_flag
    load([rootdir '/' subj '/elecs/hd_grid.mat']);
    
    if plot_dat_flag
        dat = h5read([rootdir '/DM_data/elecs.hdf5'],['/' subj]);
        dat(isnan(dat)) = 0;
        
        dat1 = dat / (max(dat) - min(dat));
        
        if strcmpi(hemi,'lh')
            offset = -5;
        else
            offset = 5;
        end
        cmap = cbrewer('seq','Reds',101);
        for i = 1:size(elecmatrix,1)
            if dat1(i) == 0
                scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
                    elecsize,'k');
            else
                scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
                    elecsize,cmap(round(dat1(i)*100)+1,:),'filled','MarkerEdgeColor','k');
            end
            hold on;
        end
    else
        if strcmpi(hemi,'lh')
            offset = -5;
        else
            offset = 5;
        end
        for i = 1:size(elecmatrix,1)
            scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
                elecsize,'y','filled','MarkerEdgeColor','k');
%             text(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%                 num2str(i),'Color','g');
            hold on;
        end
        
    end
end