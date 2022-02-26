function [nearest_chan, coords] = spmutils_channels_at_coords(SPM)
 
% Input: SPM structure with regular and custom fields
% Necessary additional SPM fields
% SPM.D:            filepath to meeg data with channel information (string)
% SPM.coords:       coordinates of interest (either in Talairach space or 
%                   voxel space)
% SPM.coord_type:   specification of coordinate type (string 'TAL' or 'VOX') 

D = spm_eeg_load(SPM.D);
coords = SPM.coords;
coord_type = SPM.coord_type;

if ~isfield(SPM,'DIM')
    DIM = [32,32,590];
elseif ~isfield(SPM,'M')
    M = [4.2500, 0, 0, -68.0000;
         0, 5.3750, 0, -100.0000;
         0, 0, 1.9531, -101.5625;
         0, 0, 0, 1.0000];
else
    DIM = SPM.xVol.DIM;
    M = SPM.xVol.M;
end

% check channels
[mod, Cind] = spm_eeg_modality_ui(D, 1, 1);
otherind = setdiff(1:nchannels(D), Cind);
if ~isempty(otherind)
    D = chantype(D, otherind, 'Other');
end
chanlabels = D.chanlabels(Cind);

% get channel locations (indices)
[chan_inds, x, y] = spm_eeg_locate_channels(D, DIM(1), Cind);

if strcmp(coord_type,'TAL')
    % convert chaninds to Talairach
    chan_inds = M * [chan_inds'; ones(2,size(chan_inds,1))];
    chan_inds = chan_inds(1:2,:)';
    pos = [chan_inds'; zeros(1,size(chan_inds,1))];
elseif strcmp(coord_type,'VOX')
    pos = [chan_inds'; zeros(1,size(chan_inds,1))];
end

% get nearest match of channel coordinates (index)
[coords,i,d] = spm_XYZreg('NearestXYZ',[coords(1); coords(2); 0],pos);

% read out index in channel list
nearest_chan = chanlabels{i};

end
