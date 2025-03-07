function failed_bug3177

% WALLTIME 00:20:00
% MEM 1gb
% DEPENDENCY ft_electroderealign mesh2edge poly2tri

%%

load(dccnpath('/project/3031000.02/external/download/tutorial/headmodel_fem/vol.mat'))


%% project electrodes
% this converts on the fly the hexaehders into a polygonal surface mesh and subsequently into a triangulated surface mesh

cfg                 = [];
cfg.elec            = ft_read_sens('standard_1020.elc');
cfg.headshape       = vol;
cfg.method          = 'project';
elec_proj           = ft_electroderealign(cfg);
