function [betas, SE, CI, nearest_chan] = spmutils_betas_at_coords(SPM)

% Input: SPM structure with regular and custom fields
% Necessary additional SPM fields
% SPM.D:            filepath to meeg data with channel information (string)
% SPM.coords:       coordinates of interest (either in Talairach space or 
%                   voxel space)
% SPM.coord_type:   specification of coordinate type (string 'TAL' or 'VOX')
% SPM.nsubs:        number of subjects
% SPM.nbetas:       number of betas

D = spm_eeg_load(SPM.D);
nsubs = SPM.nsubs;
nbetas = SPM.nbetas;
plotflag = SPM.plotflag;

% transform coords
vox_coords = SPM.coords;
if ~isfield(SPM,'coord_type')
    coord_type = 'VOX';
elseif strcmp(SPM.coord_type,'TAL')
    SPM.direction = 'TALtoVOX';
    [tal_coords, vox_coords] = spmutils_transform_coords(SPM); 
    coord_type = SPM.coord_type;
end

nearest_chan = spmutils_channels_at_coords(SPM);

CI    = 1.6449;  
eye_con = [eye(nbetas), repmat(ones(1,nsubs)*1/nsubs,nbetas,1)]';

cd(SPM.swd)
beta  = spm_data_read(SPM.Vbeta,'xyz',vox_coords);
ResMS = spm_data_read(SPM.VResMS,'xyz',vox_coords);
Bcov  = ResMS*SPM.xX.Bcov; 
betas = eye_con'*beta;
SE    = sqrt(diag(eye_con'*Bcov*eye_con));
CI    = CI*SE;

usecolor = [.4 .4 .4];
if plotflag 
    figure; hold on

    bp = bar(1:nbetas,betas,'FaceColor',usecolor,'handlevisibility','off','EdgeColor','none');
    ep = errorbar(1:nbetas,betas,CI,'Color',usecolor,'linestyle','none','handlevisibility','off');

    ep.CapSize = 1;    

    box off
    xlabel('')
    set(gca,'xtick',[])
    % set(gca,'ytick',[])
    ylabel('Beta weight')
    set(gcf,'color','white')
    ax = gca;
    ax.XColor = 'none';

    title(sprintf('Betas at %s',nearest_chan))

    
end

end