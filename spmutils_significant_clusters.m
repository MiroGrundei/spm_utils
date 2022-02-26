function [cluster, cluster_indices] = spmutils_significant_clusters(SPM)
%% Get stats from SPM contrast

% Input: SPM structure with regular and custom fields
% necessary SPM fields:
% SPM.D:            filepath to meeg data with channel and time
%                   information (string)
% SPM.usecon:       contrast of interest (is compared to the existing SPM
%                   contrasts)
% SPM.contype:      contrast type 'T' or 'F' 
% SPM.Im:           implicit masking,  empty [] or string with mask .nii
% SPM.thresDesc     'none' or 'FWE'
% SPM.u             p value threshold, e.g. 0.001;
% SPM.k:            minimum voxel extent, e.g. 0;

% optional
% SPM.coords:       coordinates of specific cluster to plot (voxel space)
% SPM.nclust:       number of cluster peaks to extract 
% SPM.plotflag:     plotting option (0=no plot; 1=all in subplots; 2=single
%                   plot for each peak

template_path = pwd;

if ~isfield(SPM,'coords')
    coords = [];
else
    coords = SPM.coords;
end

if ~isfield(SPM,'plotflag')
    plotflag = 0;
else
    plotflag = SPM.plotflag;
end

if ~isfield(SPM,'nclust')
    nclust = [];
else
    nclust = SPM.nclust;
end

contype = SPM.contype;
D = spm_eeg_load(SPM.D);
time = D.time;

% index of contrast SPM.Ic
for i = 1:numel(SPM.xCon)
    if SPM.xCon(i).c == SPM.usecon & strcmp(SPM.xCon(i).STAT,contype)
        SPM.Ic = i;
        % save for later
        Ic = i;
    end
end

% get STATS in SPM structure
[SPM,xSPM] = spm_getSPM(SPM);

% extract significant clusters
[N,Z,XYZ,A,L]  = spm_max(xSPM.Z,xSPM.XYZ);

% either plot all significant cluster peaks or only at coords of interest
if ~isempty(coords)
    % find cluster at coordinates of interest
    inclust = zeros(1,numel(L));
    for i=1:numel(L)    
        clust = L{i};
        for j=1:size(clust,2)
            if isequal(clust(:,j),coords')
                inclust(i) = 1;
            end
        end
    end
    clust_ind = find(inclust);
    if numel(clust_ind)>1
        error('More than one cluster with same coords')
    elseif numel(clust_ind) == 0
        warning('Cluster not found in stats. Showing input coordinates.')
        cluster = {coords}; 
    else
        cluster = L(clust_ind);
        Zvals = Z(A==clust_ind);
        Zcoords = XYZ(:,A==clust_ind);
    end
else
    if isempty(nclust) || nclust>length(L)
        cluster = L(1:end);
    else
        cluster = L(1:nclust);
    end
end
   
% Extract T-Maps of clusters    
   
% load topo template
template = load(fullfile(template_path,'spm_topo_template.mat'));
topo_template = template.topo_template;

% load tmap
tmap = spm_data_read(SPM.xCon(Ic).Vspm);

nc = numel(cluster);
sb_r = 2; sb_c = 2;
fgcount=1; sbcount=1;
cluster_indices = zeros(3,nc);
SPM.coord_type = 'VOX';
SPM.D = fullfile(D.path,D.fname);
for i=1:nc

    cind = cluster{i};
    time_cluster = time(cind(3,:));
    
    temp = [];
    for j=1:size(cind,2)
        temp = [temp,tmap(cind(1,j),cind(2,j),cind(3,j))];
    end
    [m, max_ind] = max(temp);
    
    cluster_indices(:,i) = cind(:,max_ind);
    
    SPM.coords = cluster_indices(:,i);
    nearest_chan = spmutils_channels_at_coords(SPM);

    if time(cind(3)) > 0 % only plot if nonnegative time
       
    if plotflag
        
        plotmap = tmap(:,:,cind(3,max_ind));       
        Y = fliplr(rot90(rot90(rot90(plotmap+topo_template))));
        
        if plotflag == 1
            figure(fgcount); subplot(sb_r,sb_c,sbcount);
            sbcount=sbcount+1;
            if sbcount==sb_r*sb_c +1
                fgcount=fgcount+1;
                sbcount=1;
            end
        else
            figure(i)
        end
        surface(Y,'edgecolor','none');
        colormap('jet')
        h = colorbar;
        axis square
        axis off
        box off
        ylabel(h, sprintf('%s value',contype))
        set(gcf,'color','w');

        title(sprintf('Peak at %.0fms at %s (Cluster: %.0f-%.0fms)', ...
                        1000*time(cind(3,max_ind)), nearest_chan, 1000*time_cluster(1),  1000*time_cluster(end)))
        
    end
    
    end % only plot if non negative time
end

% Show peak times
fprintf('Cluster peaks at: ')
for t=1:size(cluster_indices,2)
    fprintf('%.0fms ' , 1000*D.time(cluster_indices(3,t)))
end
fprintf('\n')

end
    
    
