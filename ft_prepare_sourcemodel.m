function [sourcemodel, cfg] = ft_prepare_sourcemodel(cfg)

% FT_PREPARE_SOURCEMODEL constructs a source model, for example a 3D grid or a
% cortical sheet. The source model that can be used for source reconstruction,
% beamformer scanning, linear estimation and MEG interpolation.
%
% Use as
%   sourcemodel = ft_prepare_sourcemodel(cfg)
% where the details of the configuration structure determine how the source
% model will be constructed.
%
% The different approaches for constructing a source model are
%   cfg.method = 'basedongrid'        regular 3D grid with explicit specification
%                'basedonresolution'  regular 3D grid with specification of the resolution
%                'basedonpos'         place dipoles at the predefined positions
%                'basedonmri'         regular 3D grid, based on segmented MRI, restricted to gray matter
%                'basedonmni'         regular 3D grid, based on a warped template grid, based on the MNI brain
%                'basedoncortex'      cortical sheet from external software such as Caret or FreeSurfer, can also be two separate hemispheres
%                'basedonshape'       surface mesh based on inward shifted head surface from an external file
%                'basedonvol'         surface mesh based on inward shifted brain surface from volume conductor
%                'basedonfile'        the sourcemodel should be read from file
%                'basedoncentroids'   irregular 3D grid based on volumetric mesh
% The default method is determined automatically based on the configuration options
% that you specify.
%
% BASEDONGRID - uses an explicitly specified grid, according to the following
% configuration options:
%   cfg.xgrid         = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.ygrid         = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.zgrid         = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%
% BASEDONRESOLUTION - uses an grid with the desired resolution, according
% to the following configuration options:
%   cfg.resolution    = number (e.g. 1 cm) for automatic grid generation
%
% BASEDONPOS - places sources on positions that you explicitly specify, according to
% the following configuration options:
%   cfg.sourcemodel.pos       = N*3 matrix with position of each source
%   cfg.sourcemodel.inside    = N*1 vector with boolean value whether position is inside brain (optional)
%   cfg.sourcemodel.dim       = [Nx Ny Nz] vector with dimensions in case of 3D grid (optional)
% The following fields (from FT_PRERARE_LEADFIELD or FT_SOURCEANALYSIS) are
% not used in this function, but will be copied along to the output:
%   cfg.sourcemodel.leadfield = cell-array
%   cfg.sourcemodel.filter    = cell-array
%   cfg.sourcemodel.subspace
%   cfg.sourcemodel.lbex
%
% BASEDONMNI - uses source positions from a template sourcemodel that is inversely
% warped from MNI coordinates to the individual subjects MRI. It uses the following
% configuration options:
%   cfg.mri             = structure with the anatomical MRI, or the filename of the MRI, see FT_READ_MRI
%   cfg.nonlinear       = 'no' (or 'yes'), use non-linear normalization
%   cfg.resolution      = number (e.g. 6) of the resolution of the template MNI grid, defined in mm
%   cfg.template        = structure with the template sourcemodel, or the filename of a template sourcemodel (defined in MNI space)
%   cfg.templatemri     = string, filename of the MNI template (default = 'T1.mnc' for SPM2 or 'T1.nii' for SPM8 and SPM12)
%   cfg.spmversion      = string, 'spm2', 'spm8', 'spm12' (default = 'spm12')
%   cfg.spmmethod       = string, 'old', 'new' or 'mars', see FT_VOLUMENORMALISE
%   cfg.nonlinear       = string, 'yes' or 'no', see FT_VOLUMENORMALISE
% Either cfg.resolution or cfg.template needs to be defined; if both are defined, cfg.template prevails.
%
% BASEDONMRI - makes a segmentation of the individual anatomical MRI and places
% sources in the grey matter. It uses the following configuration options:
%   cfg.mri             = can be filename, MRI structure or segmented MRI structure
%   cfg.threshold       = 0.1, relative to the maximum value in the segmentation
%   cfg.smooth          = 5, smoothing in voxels
%
% BASEDONCORTEX - places sources on the vertices of a cortical surface description
%   cfg.headshape       = string, should be a *.fif file
%
% BASEDONCENTROIDS - places sources on the centroids of a volumetric mesh
%   cfg.headmodel       = tetrahedral or hexahedral mesh
%   cfg.headmodel.type  = 'simbio'
%
% Other configuration options include
%   cfg.unit            = string, can be 'mm', 'cm', 'm' (default is automatic, based on the input data)
%   cfg.tight           = 'yes' or 'no' (default is automatic)
%   cfg.inwardshift     = number, amount to shift the innermost surface of the headmodel inward when determining
%                         whether sources are inside or outside the source compartment (default = 0)
%   cfg.moveinward      = number, amount to move sources inward to ensure a certain minimal distance to the innermost
%                         surface of the headmodel (default = 0)
%   cfg.movetocentroids = 'yes' or 'no', move the dipoles to the centroids of the hexahedral
%                         or tetrahedral mesh (default = 'no')
%   cfg.spherify        = 'yes' or 'no', scale the source model so that it fits inside a sperical
%                         volume conduction model (default = 'no')
%   cfg.symmetry        = 'x', 'y' or 'z' symmetry for two dipoles, can be empty (default = [])
%   cfg.headshape       = a filename for the headshape, a structure containing a single surface,
%                         or a Nx3 matrix with headshape surface points (default = [])
%   cfg.spmversion      = string, 'spm2', 'spm8', 'spm12' (default = 'spm12')
%
% The EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad          = structure with gradiometer definition or filename, see FT_READ_SENS
%
% The headmodel or volume conduction model can be specified as
%   cfg.headmodel     = structure with volume conduction model or filename, see FT_PREPARE_HEADMODEL
%
% The cfg.inwardshift option can be used for 3D grids to specify a positive (inward)
% or negative (outward) number to shift the innermost surface of the headmodel
% (usually the skull) when determining whether sources are to be flagged as inside or
% outside the source compartment. Only sources flagged as inside will be considered
% for subsequent source reconstructions. An ourward shift can be useful for a
% spherical or singleshell MEG headmodel. For a source model based on a cortical
% sheet in general you want all sources to be considered inside. For a BEM headmodel
% (EEG or MEG), there should never be any sources outside the actual source
% compartment.
%
% The cfg.moveinward option can be used for a source model based on a cortical sheet
% to push the sources inward a little bit to ensure sufficient distance to the
% innermost surface of a BEM headmodel (EEG or MEG).
%
% See also FT_PREPARE_LEADFIELD, FT_PREPARE_HEADMODEL, FT_SOURCEANALYSIS,
% FT_DIPOLEFITTING, FT_MEGREALIGN

% Copyright (C) 2004-2024, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance headmodel sens

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', 'mriunits');
cfg = ft_checkconfig(cfg, 'forbidden', 'sourceunits');
cfg = ft_checkconfig(cfg, 'renamed', {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'vol',     'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'grid',    'sourcemodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed', {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed', {'optofile', 'opto'});
cfg = ft_checkconfig(cfg, 'renamedval', {'unit', 'auto', []});

cfg = ft_checkconfig(cfg, 'renamed', {'tightgrid', 'tight'});  % this is moved to cfg.sourcemodel.tight by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'renamed', {'sourceunits', 'unit'}); % this is moved to cfg.unit by the subsequent createsubcfg

% put the low-level options pertaining to the sourcemodel in their own field
cfg = ft_checkconfig(cfg, 'createsubcfg', {'sourcemodel'});
% move some fields from cfg.sourcemodel back to the top-level configuration
cfg = ft_checkconfig(cfg, 'createtopcfg', {'sourcemodel'});

% set the defaults
cfg.method            = ft_getopt(cfg, 'method'); % the default is to do automatic detection further down
cfg.headshape         = ft_getopt(cfg, 'headshape');
cfg.mri               = ft_getopt(cfg, 'mri');
cfg.headmodel         = ft_getopt(cfg, 'headmodel');
cfg.sourcemodel       = ft_getopt(cfg, 'sourcemodel');
cfg.unit              = ft_getopt(cfg, 'unit');
cfg.symmetry          = ft_getopt(cfg, 'symmetry');
cfg.spmversion        = ft_getopt(cfg, 'spmversion', 'spm12');
cfg.spherify          = ft_getopt(cfg, 'spherify', 'no');
cfg.movetocentroids   = ft_getopt(cfg, 'movetocentroids', 'no');
cfg.moveinward        = ft_getopt(cfg, 'moveinward'); % the default is automatic and depends on a triangulation being present
cfg.checkinside       = ft_getopt(cfg, 'checkinside', 'no'); % default is 'no' since this is a relatively slow procedure. It is also not always required, for example with MEG singlesphere, singleshell, localspheres.
cfg.feedback          = ft_getopt(cfg, 'feedback', 'text');

% this option was deprecated on 12 Aug 2020
if isfield(cfg, 'warpmni')
  % prior to the introduction of cfg.method we used cfg.warpmni to separate between
  % basedonmni and basedonmri, which both require cfg.mri to be present
  if isfield(cfg, 'mri') && istrue(cfg.warpmni)
    ft_warning('please specify cfg.method=''basedonmni'' instead of cfg.warpmni=''yes''');
    cfg.method = 'basedonmni';
  elseif isfield(cfg, 'mri') && ~istrue(cfg.warpmni)
    ft_warning('please specify cfg.method=''basedonmri'' instead of cfg.warpmni=''no''');
    cfg.method = 'basedonmri';
  end
  cfg = rmfield(cfg, 'warpmni');
end

% this code expects the inside to be represented as a logical array
cfg = ft_checkconfig(cfg, 'inside2logical', 'yes');
if isfield(cfg, 'sourcemodel')
  cfg.sourcemodel = ft_checkconfig(cfg.sourcemodel, 'renamed',  {'pnt' 'pos'});
end
if isfield(cfg, 'template')
  cfg.template = ft_checkconfig(cfg.template, 'renamed',  {'pnt' 'pos'});
end

if isfield(cfg, 'resolution') && isfield(cfg, 'xgrid') && ~ischar(cfg.xgrid)
  ft_error('You cannot specify cfg.resolution and an explicit cfg.xgrid simultaneously');
end
if isfield(cfg, 'resolution') && isfield(cfg, 'ygrid') && ~ischar(cfg.ygrid)
  ft_error('You cannot specify cfg.resolution and an explicit cfg.ygrid simultaneously');
end
if isfield(cfg, 'resolution') && isfield(cfg, 'zgrid') && ~ischar(cfg.zgrid)
  ft_error('You cannot specify cfg.resolution and an explicit cfg.zgrid simultaneously');
end

% the source model can be constructed in a number of ways
if isempty(cfg.method)
  if isfield(cfg, 'sourcemodel') && ischar(cfg.sourcemodel)
    cfg.method = 'basedonfile';
  elseif isfield(cfg, 'xgrid') && ~ischar(cfg.xgrid)
    cfg.method = 'basedongrid'; % regular 3D grid with explicit specification
  elseif isfield(cfg.sourcemodel, 'pos')
    cfg.method = 'basedonpos'; % using user-supplied positions, which can be regular or irregular
  elseif ~isempty(cfg.headshape)
    cfg.method = 'basedonshape'; % surface mesh based on inward shifted head surface from external file
  elseif ~isempty(cfg.mri)
    cfg.method = 'basedonmri'; % regular 3D grid, based on segmented MRI, restricted to gray matter
  elseif isfield(cfg, 'headshape') && (iscell(cfg.headshape) || any(ft_filetype(cfg.headshape, {'neuromag_fif', 'freesurfer_triangle_binary', 'caret_surf', 'gifti'})))
    cfg.method = 'basedoncortex'; % cortical sheet from external software such as Caret or FreeSurfer, can also be two separate hemispheres
  elseif isfield(cfg, 'resolution')
    cfg.method = 'basedonresolution'; % regular 3D grid with specification of the resolution
  elseif ~isempty(cfg.headmodel)
    cfg.method = 'basedonvol'; % surface mesh based on inward shifted brain surface from volume conductor
  else
    ft_error('incorrect cfg specification for constructing a sourcemodel');
  end
else
  cfg = ft_checkconfig(cfg, 'allowedval', {'method', 'basedongrid', 'basedonpos', 'basedonshape', ...
    'basedonmri', 'basedonmni', 'basedoncortex', 'basedonresolution', 'basedonvol', 'basedonfile','basedoncentroids'});
end

% these are mutually exclusive, but printing all requested methods here
% facilitates debugging of weird configs. Also specify the defaults here to
% keep the overview
switch cfg.method
  case 'basedonfile'
    fprintf('reading sourcemodel from file\n');
    cfg.tight       = ft_getopt(cfg, 'tight',   'no');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 0); % in this case for inside detection

  case 'basedongrid'
    fprintf('creating sourcemodel based on user-specified 3D grid\n');
    cfg.xgrid       = ft_getopt(cfg, 'xgrid', 'auto');
    cfg.ygrid       = ft_getopt(cfg, 'ygrid', 'auto');
    cfg.zgrid       = ft_getopt(cfg, 'zgrid', 'auto');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 0); % in this case for inside detection
    cfg.tight       = ft_getopt(cfg, 'tight',   'yes');
    cfg.resolution  = [];

  case 'basedonresolution'
    fprintf('creating sourcemodel based on automatic 3D grid with the specified resolution\n');
    cfg.xgrid       = ft_getopt(cfg, 'xgrid', 'auto');
    cfg.ygrid       = ft_getopt(cfg, 'ygrid', 'auto');
    cfg.zgrid       = ft_getopt(cfg, 'zgrid', 'auto');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 0); % in this case for inside detection
    if istrue(cfg.movetocentroids)
      cfg.tight     = ft_getopt(cfg, 'tight',   'no');
    else
      cfg.tight     = ft_getopt(cfg, 'tight',   'yes');
    end

  case 'basedonpos'
    fprintf('creating sourcemodel based on user specified dipole positions\n');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 0); % in this case for inside detection
    cfg.tight       = ft_getopt(cfg, 'tight',    'no');

  case 'basedonshape'
    fprintf('creating sourcemodel based on inward-shifted head shape\n');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift',  0); % in this case for inside detection
    cfg.spheremesh  = ft_getopt(cfg, 'spheremesh', 642);
    cfg.tight       = ft_getopt(cfg, 'tight',    'yes');

  case 'basedoncortex'
    fprintf('creating sourcemodel based on cortical sheet from an external file\n');
    cfg.tight       = ft_getopt(cfg, 'tight',     'no');

  case 'basedonmri'
    fprintf('creating sourcemodel based on an anatomical volume\n');
    cfg.threshold   = ft_getopt(cfg, 'threshold', 0.1); % relative
    cfg.smooth      = ft_getopt(cfg, 'smooth',      5); % in voxels
    cfg.tight       = ft_getopt(cfg, 'tight',   'yes');

  case 'basedonvol'
    fprintf('creating sourcemodel based on inward-shifted brain surface from volume conductor model\n');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift',   0); % in this case for inside detection
    cfg.spheremesh  = ft_getopt(cfg, 'spheremesh',  642);
    cfg.tight       = ft_getopt(cfg, 'tight',      'no');

  case 'basedonmni'
    cfg.tight       = ft_getopt(cfg.sourcemodel, 'tight',       'no');
    cfg.nonlinear   = ft_getopt(cfg.sourcemodel, 'nonlinear',   'no');

  case 'basedoncentroids'
    fprintf('creating sourcemodel based on volumetric mesh centroids\n');
    cfg.tight       = ft_getopt(cfg.sourcemodel, 'tight',       'no');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 0); % in this case for inside detection
end

if (isfield(cfg, 'smooth') && ~strcmp(cfg.smooth, 'no')) || strcmp(cfg.method, 'basedonmni')
  % check that the preferred SPM version is on the path
  ft_hastoolbox(cfg.spmversion, 1);
end

% start with an empty structure
sourcemodel = [];

% get the volume conduction model
if ~isempty(cfg.headmodel) && ischar(cfg.headmodel)
  headmodel = ft_read_headmodel(cfg.headmodel);
else
  % ensure that the volume conduction model is up-to-date
  headmodel = ft_datatype_headmodel(cfg.headmodel);
end

% get the headshape, this can also be a cortical sheet, or a set of left and right cortical sheets
if ~isempty(cfg.headshape) && ischar(cfg.headshape)
  headshape = ft_read_headshape(cfg.headshape);
else
  headshape = cfg.headshape;
end

% get the anatomical MRI or segmentation
if ~isempty(cfg.mri) && ischar(cfg.mri)
  mri = ft_read_mri(cfg.mri);
else
  mri = cfg.mri;
end

% get the gradiometer or electrode definition
try
  sens = ft_fetch_sens(cfg);
catch
  sens = [];
end

if strcmp(cfg.method, 'basedonfile')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the source model from a MATLAB file
  % this needs to be done prior to determining the default units
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cfg.sourcemodel = loadvar(cfg.sourcemodel, 'sourcemodel');
  sourcemodel = cfg.sourcemodel;
end

if isempty(cfg.unit)
  if isfield(cfg.sourcemodel, 'unit') && ~isempty(cfg.sourcemodel.unit)
    % take the existing source model units
    cfg.unit = cfg.sourcemodel.unit;
  elseif isfield(cfg.sourcemodel, 'pos') && size(cfg.sourcemodel.pos,1)>10
    % estimate the units based on the existing source positions
    cfg.sourcemodel = ft_determine_units(cfg.sourcemodel);
    cfg.unit = cfg.sourcemodel.unit;
  elseif strcmp(cfg.method, 'basedonmni') && ~isempty(cfg.mri.unit)
    % take the existing MRI units
    cfg.unit = cfg.mri.unit;
  elseif strcmp(cfg.method, 'basedonmri') && ~isempty(cfg.mri.unit)
    % take the existing MRI units
    cfg.unit = cfg.mri.unit;
  elseif ~isempty(sens)
    % take the units from the gradiometer or electrode array
    cfg.unit = sens.unit;
  elseif ~isempty(headmodel)
    % take the units from the volume conduction model
    cfg.unit = headmodel.unit;
  else
    ft_warning('assuming "cm" as default units for source model');
    cfg.unit = 'cm';
  end
  ft_warning('assuming that the sourcemodel units are in %s', cfg.unit);
end

% convert the volume conduction model to the desired units for the source model
if ~isempty(headmodel)
  headmodel = ft_convert_units(headmodel, cfg.unit);
end

% convert the headshape to the desired units for the source model
if ~isempty(headshape)
  headshape = ft_convert_units(headshape, cfg.unit);
end

% convert the sensor array to the desired units for the source model
if ~isempty(sens)
  sens = ft_convert_units(sens, cfg.unit);
end

if ~isempty(mri)
  % convert the mri to the desired units for the source model
  mri = ft_convert_units(mri, cfg.unit);
end

switch cfg.method
  case 'basedongrid'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % a detailed xgrid/ygrid/zgrid has been specified, the other details
    % still need to be determined
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ischar(cfg.xgrid) || ischar(cfg.ygrid) || ischar(cfg.zgrid)
      % th especification 'auto' is not allowed here
      ft_error('you must specify explicit values for xgrid/ygrid/zgrid')
    end
    sourcemodel.dim   = [length(cfg.xgrid) length(cfg.ygrid) length(cfg.zgrid)];
    [X, Y, Z]  = ndgrid(cfg.xgrid, cfg.ygrid, cfg.zgrid);
    sourcemodel.pos   = [X(:) Y(:) Z(:)];
    sourcemodel.unit  = cfg.unit;

  case 'basedonresolution'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct a regular 3D grid that spans a box encompassing the complete brain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('creating 3D grid with %g %s resolution\n', cfg.resolution, cfg.unit);

    if ~isempty(sens) && isfield(sens, 'chanpos')
      % determine the bounding box of the sensor array
      minpos = min(sens.chanpos,[],1);
      maxpos = max(sens.chanpos,[],1);
    elseif ~isempty(headmodel)
      % determine the bounding box of the volume conduction model
      if isfield(headmodel, 'bnd') && isfield(headmodel.bnd, 'pos')
        pos = cat(1, headmodel.bnd(:).pos);
      elseif isfield(headmodel, 'bnd') && isfield(headmodel.bnd, 'pnt')
        pos = cat(1, headmodel.bnd(:).pnt);
      elseif isfield(headmodel, 'pos')
        pos = headmodel.pos;
      elseif ft_headmodeltype(headmodel, 'localspheres')
        pos = headsurface(headmodel, sens);
      elseif ft_headmodeltype(headmodel, 'singlesphere')
        pos = [
          headmodel.o - headmodel.r
          headmodel.o + headmodel.r
          ];
      elseif ft_headmodeltype(headmodel, 'concentricspheres')
        pos = [
          headmodel.o - max(headmodel.r)
          headmodel.o + max(headmodel.r)
          ];
      end
      minpos = min(pos,[],1);
      maxpos = max(pos,[],1);
      % add a few percent on either side
      minpos(minpos<0) = minpos(minpos<0).*1.08;
      maxpos(maxpos>0) = maxpos(maxpos>0).*1.08;
      minpos(minpos>0) = minpos(minpos>0).*0.92;
      maxpos(maxpos<0) = maxpos(maxpos<0).*0.92;
    elseif ~isempty(headshape)
      % determine the bounding box of the headshape
      minpos = min(headshape.pos,[],1);
      maxpos = max(headshape.pos,[],1);
      % add a few percent on either side
      minpos(minpos<0) = minpos(minpos<0).*1.08;
      maxpos(maxpos>0) = maxpos(maxpos>0).*1.08;
      minpos(minpos>0) = minpos(minpos>0).*0.92;
      maxpos(maxpos<0) = maxpos(maxpos<0).*0.92;
    else
      ft_error('creating an automatic 3D grid requires either the sensor positions, a headmodel, or a headshape to estimate the extent');
    end

    if isempty(cfg.symmetry)
      % round the limits such that [0 0 0] will be on the grid
      minpos = floor(minpos/cfg.resolution)*cfg.resolution;
      maxpos = ceil(maxpos/cfg.resolution)*cfg.resolution;
    else
      % round the limits such that the grid will be symmetric around [0 0 0]
      minpos = floor((minpos+cfg.resolution/2)/cfg.resolution)*cfg.resolution - cfg.resolution/2;
      maxpos = ceil((maxpos+cfg.resolution/2)/cfg.resolution)*cfg.resolution - cfg.resolution/2;
    end

    if ischar(cfg.xgrid) && strcmp(cfg.xgrid, 'auto')
      xgrid = minpos(1):cfg.resolution:maxpos(1);
    else
      xgrid = cfg.xgrid(1):cfg.resolution:cfg.xgrid(end);
    end
    if ischar(cfg.ygrid) && strcmp(cfg.ygrid, 'auto')
      ygrid = minpos(2):cfg.resolution:maxpos(2);
    else
      ygrid = cfg.ygrid(1):cfg.resolution:cfg.ygrid(end);
    end
    if ischar(cfg.zgrid) && strcmp(cfg.zgrid, 'auto')
      zgrid = minpos(3):cfg.resolution:maxpos(3);
    else
      zgrid = cfg.zgrid(1):cfg.resolution:cfg.zgrid(end);
    end
    sourcemodel.dim   = [length(xgrid) length(ygrid) length(zgrid)];
    [X, Y, Z]  = ndgrid(xgrid, ygrid, zgrid);
    sourcemodel.pos   = [X(:) Y(:) Z(:)];
    sourcemodel.unit  = cfg.unit;
    fprintf('initial 3D grid dimensions are [%d %d %d]\n', sourcemodel.dim(1), sourcemodel.dim(2), sourcemodel.dim(3));

  case 'basedonpos'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % source positions are already specified in the configuration, reuse as much of the
    % prespecified model as possible (but only known objects)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sourcemodel = keepfields(cfg.sourcemodel, {'pos', 'tri', 'dim', 'transform', 'unit', 'coordsys', 'xgrid', 'ygrid', 'zgrid', 'mom', 'inside', 'lbex', 'subspace', 'leadfield', 'leadfielddimord', 'filter', 'filterdimord', 'label'});

  case 'basedonmri'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct a grid based on the segmented MRI that is provided in the
    % configuration, only voxels in gray matter will be used
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isfield(cfg, 'resolution')
      switch cfg.unit
        case 'mm'
          cfg.resolution = 10;
        case 'cm'
          cfg.resolution = 1;
        case 'dm'
          cfg.resolution = 0.1;
        case 'm'
          cfg.resolution = 0.01;
      end
    end

    issegmentation = false;
    if isfield(mri, 'gray')
      % this is based on tissue probability maps, being the original implementation here.
      dat = double(mri.gray);

      % apply a smoothing of a certain amount of voxels
      if ~strcmp(cfg.smooth, 'no')
        dat = volumesmooth(dat, cfg.smooth, 'MRI gray matter');
      end

    elseif isfield(mri, 'anatomy')
      % this could be an anatomical MRI but also a segmentation or tpm stored as a
      % NIFTI file, reading it from disk always leads to the field 'anatomy'.
      dat = double(mri.anatomy);

      % apply a smoothing of a certain amount of voxels
      if ~strcmp(cfg.smooth, 'no')
        dat = volumesmooth(dat, cfg.smooth, 'anatomy');
      end

    elseif ft_datatype(mri, 'segmentation')
      % this is a proper segmentation, where a set of boolean masks is in the
      % input, or and indexed volume, along with labels. FIXME for now still
      % only works for boolean volumes.
      issegmentation = true;
      fn = booleanfields(mri);
      if isempty(fn)
        % convert indexed segmentation into probabilistic
        mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic');
        fn  = booleanfields(mri);
      end

      dat = false(mri.dim);
      for i=1:numel(fn)
        if ~strcmp(cfg.smooth, 'no')
          mri.(fn{i}) = volumesmooth(double(mri.(fn{i})), cfg.smooth, fn{i}) > cfg.threshold;
        end
        dat = dat | mri.(fn{i});
      end
      dat = double(dat);
    else
      ft_error('cannot determine the format of the segmentation in cfg.mri');
    end

    % determine for each voxel whether it belongs to the grey matter
    fprintf('thresholding MRI data at a relative value of %f\n', cfg.threshold);
    head = dat./max(dat(:)) > cfg.threshold;

    ind                 = find(head(:));
    fprintf('%d from %d voxels in the segmentation are marked as ''inside'' (%.0f%%)\n', length(ind), numel(head), 100*length(ind)/numel(head));
    [X,Y,Z]             = ndgrid(1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3));             % create the grid in MRI-coordinates
    posmri              = [X(ind) Y(ind) Z(ind)];                                       % take only the inside voxels
    poshead             = ft_warp_apply(mri.transform, posmri);                         % transform to head coordinates
    resolution          = cfg.resolution;                                               % source and mri are expressed in the same units
    
    % determine the size of a bounding box, round it off to the nearest mm
    scale = ft_scalingfactor(cfg.unit, 'mm');
    xmin = floor(min(scale*poshead(:,1)))/scale;
    xmax = ceil (max(scale*poshead(:,1)))/scale;
    ymin = floor(min(scale*poshead(:,2)))/scale;
    ymax = ceil (max(scale*poshead(:,2)))/scale;
    zmin = floor(min(scale*poshead(:,3)))/scale;
    zmax = ceil (max(scale*poshead(:,3)))/scale;

    xgrid               = xmin:resolution:xmax;  % create the grid in head-coordinates
    ygrid               = ymin:resolution:ymax;  % with consistent x,y,z definitions
    zgrid               = zmin:resolution:zmax;
    [X,Y,Z]             = ndgrid(xgrid,ygrid,zgrid);
    pos2head            = [X(:) Y(:) Z(:)];
    pos2mri             = ft_warp_apply(inv(mri.transform), pos2head);                  % transform to MRI voxel coordinates
    pos2mri             = round(pos2mri);
    inside              = getinside(pos2mri, head);                                     % use helper subfunction

    sourcemodel.pos     = pos2head;
    sourcemodel.dim     = [length(xgrid) length(ygrid) length(zgrid)];
    sourcemodel.inside  = inside(:);
    sourcemodel.unit    = cfg.unit;

    if issegmentation
      % pass on the segmentation information on the grid points, the
      % individual masks have been smoothed above
      fn = booleanfields(mri);
      for i=1:numel(fn)
        sourcemodel.(fn{i}) = getinside(pos2mri, mri.(fn{i}));
      end
      % convert back is not in general possible because the masks can be
      % overlapping due to smoothing
      % sourcemodel = ft_datatype_segmentation(sourcemodel, 'segmentationstyle', segstyle);
    end
    fprintf('the full grid contains %d points\n', numel(sourcemodel.inside));
    fprintf('%d grid points are marked as inside the brain\n',  sum(sourcemodel.inside));

  case 'basedoncortex'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the headshape from a *.fif file that was created using Freesurfer and MNE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % return the vertices and triangles from the cortical sheet
    sourcemodel.pos  = headshape.pos;
    sourcemodel.tri  = headshape.tri;
    sourcemodel.unit = headshape.unit;

  case 'basedonshape'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use the headshape  to make a superficial dipole layer (e.g.
    % for megrealign). Assume that all points are inside the volume.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the surface describing the head shape
    [headshape.pos, headshape.tri] = headsurface([], [], 'headshape', cfg.headshape);

    % ensure that the headshape is in the same units as the source model
    headshape = ft_convert_units(headshape, cfg.unit);

    % note that cfg.inwardshift should be expressed in the units consistent with the data
    sourcemodel.pos     = headsurface([], [], 'headshape', headshape, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
    sourcemodel.tri     = headshape.tri;
    sourcemodel.unit    = headshape.unit;
    sourcemodel.inside  = true(size(sourcemodel.pos,1),1);

  case 'basedonvol'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use the volume conduction model to make a superficial dipole layer (e.g.
    % for megrealign). Assume that all points are inside the volume.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % please note that cfg.inwardshift should be expressed in the units consistent with cfg.unit
    sourcemodel.pos     = headsurface(headmodel, sens, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
    sourcemodel.unit    = cfg.unit;
    sourcemodel.inside  = true(size(sourcemodel.pos,1),1);

  case 'basedonmni'
    if ~isfield(cfg, 'template') && ~isfield(cfg, 'resolution')
      ft_error('you either need to specify the filename of a template grid in cfg.template, or a resolution in cfg.resolution');
    elseif isfield(cfg, 'template')
      % let the template filename prevail
      fname = cfg.template;
    elseif isfield(cfg, 'resolution')
      % use one of the templates that are in Fieldtrip, this requires a resolution
      if isequal(cfg.resolution, round(cfg.resolution))
        fname = sprintf('standard_sourcemodel3d%dmm.mat', cfg.resolution);
      else
        fname = sprintf('standard_sourcemodel3d%dpoint%dmm.mat', floor(cfg.resolution), round(10*(cfg.resolution-floor(cfg.resolution))));
      end
      if ~exist(fname, 'file')
        ft_error('the MNI template grid based on the specified resolution does not exist');
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check whether the mni template grid exists for the specified resolution
    % if not create it: FIXME (this needs to be done still)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get the template grid
    if ischar(fname)
      mnigrid = loadvar(fname, 'sourcemodel');
    else
      mnigrid = cfg.template;
    end

    % ensure these to have units in mm, the conversion of the source model into the desired units is done further down
    mri     = ft_convert_units(mri,     'mm');
    mnigrid = ft_convert_units(mnigrid, 'mm');

    % ensure that it is specified with logical inside
    mnigrid = fixinside(mnigrid);

    % spatial normalisation of the MRI to the template
    tmpcfg = keepfields(cfg, {'spmversion', 'spmmethod', 'nonlinear'});
    if isfield(cfg, 'templatemri')
      % this option is called differently for the two functions
      tmpcfg.template = cfg.templatemri;
    end
    normalise = ft_volumenormalise(tmpcfg, mri);

    % the normalisation from original subject head coordinates to MNI consists of an initial rigid body transformation, followed by a more precise (non)linear transformation
    % the reverse transformations are used to get from MNI to the original subject head coordinates
    % first apply the inverse of the nonlinear transformation, followed by the inverse of the initial rigid body transformation
    sourcemodel.pos = ft_warp_apply(inv(normalise.initial), ft_warp_apply(normalise.params, mnigrid.pos, 'sn2individual'), 'homogeneous');

    % copy some of the fields over from the input arguments
    sourcemodel = copyfields(mri,       sourcemodel, {'unit', 'coordsys'});
    sourcemodel = copyfields(mnigrid,   sourcemodel, {'dim', 'tri', 'inside'});
    sourcemodel = copyfields(normalise, sourcemodel, {'params', 'initial'});
    if ft_datatype(mnigrid, 'parcellation')
      % copy the boolean fields over from the template MNI grid
      sourcemodel = copyfields(mnigrid, sourcemodel, booleanfields(mnigrid));
    end

  case 'basedoncentroids'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the centroids of each volume element of a FEM mesh
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this will also copy some fields from the headmodel, such as tissue, tissuelabel, unit, coordsys
    sourcemodel = compute_centroids(headmodel);
end

if isfield(sourcemodel, 'unit')
  % in most cases the source model will already be in the desired units, but for
  % "basedonmni" it will be in 'mm' since in the spatial normalization the MRI and
  % the MNI template are converted to 'mm'
  sourcemodel = ft_convert_units(sourcemodel, cfg.unit);
else
  % the units were specified by the user or determined automatically, assign them to the source model
  sourcemodel.unit = cfg.unit;
end

% do some sanity checks
if isfield(sourcemodel, 'filter')
  assert(numel(sourcemodel.filter) == size(sourcemodel.pos, 1), 'the number of precomputed filters does not match number of source positions');
end
if isfield(sourcemodel, 'leadfield')
  assert(numel(sourcemodel.leadfield) == size(sourcemodel.pos, 1), 'the number of precomputed leadfields does not match number of source positions');
end

if strcmp(cfg.spherify, 'yes')
  if ~ft_headmodeltype(headmodel, 'singlesphere') && ~ft_headmodeltype(headmodel, 'concentricspheres')
    ft_error('this only works for spherical volume conduction models');
  end
  % deform the cortex so that it fits in a unit sphere
  pos = mesh_spherify(sourcemodel.pos, [], 'shift', 'range');
  % scale it to the radius of the innermost sphere, make it a tiny bit smaller to
  % ensure that the support point with the exact radius 1 is still inside the sphere
  pos = pos*min(headmodel.r)*0.999;
  pos(:,1) = pos(:,1) + headmodel.o(1);
  pos(:,2) = pos(:,2) + headmodel.o(2);
  pos(:,3) = pos(:,3) + headmodel.o(3);
  sourcemodel.pos = pos;
end % if spherify

if ~isempty(cfg.moveinward)
  if ~ismember(cfg.method, {'basedonshape', 'basedoncortex', 'basedonvol', 'basedonfile'})
    ft_warning('cfg.moveinward is designed to work with surface based sourcemodels, not with 3D grid sourcemodels.')
  end
  % construct a triangulated boundary of the source compartment
  [pos1, tri1] = headsurface(headmodel, [], 'inwardshift', cfg.moveinward, 'surface', 'brain');
  inside = surface_inside(sourcemodel.pos, pos1, tri1);
  if ~all(inside)
    pos2 = sourcemodel.pos(~inside,:);
    [dum, pos3] = project_elec(pos2, pos1, tri1);
    sourcemodel.pos(~inside,:) = pos3;
  end
  if cfg.moveinward>cfg.inwardshift
    sourcemodel.inside  = true(size(sourcemodel.pos,1),1);
  end
end % if moveinward

if strcmp(cfg.movetocentroids, 'yes')
  % compute centroids of the tetrahedral or hexahedral mesh
  centroids = compute_centroids(headmodel);

  % move the dipole positions in the sourcemodel to the closest centroid
  indx = knnsearch(centroids.pos, sourcemodel.pos);
  sourcemodel.pos = centroids.pos(indx,:);

  if isfield(centroids, 'tissue') && isfield(centroids, 'tissuelabel')
    % remember the tissue type around each dipole
    sourcemodel.tissue = centroids.tissue(indx);
    sourcemodel.tissuelabel = centroids.tissuelabel;
  end

  % these fields can be copied over from the headmodel
  sourcemodel = copyfields(centroids, sourcemodel, {'coordsys', 'unit'});

  % eliminate duplicate positions, this applies for example if cfg.resolution is smaller than the mesh resolution
  [sourcemodel.pos, indx] = unique(sourcemodel.pos, 'rows', 'stable');
  % also eliminate the duplicates in the tissue
  if isfield(sourcemodel, 'tissue') && isfield(sourcemodel, 'tissuelabel')
    sourcemodel.tissue = sourcemodel.tissue(indx);
  end

  % the shifted positions are not on a regular 3D grid any more, hence dim does not apply
  sourcemodel = removefields(sourcemodel, {'dim'});
end % if movetocentroids

if isfield(sourcemodel, 'inside') && isfield(cfg, 'inwardshift') && isfield(cfg, 'template')
  % warn about inwardshift not having an effect as inside is already specified as well
  % warning should only be issued for templates, inwardshift can also be present for surface meshes
  ft_warning('Inside dipole locations already determined by a template, cfg.inwardshift has no effect.')
end

% determine the dipole locations that are inside the source compartment of the
% volume conduction model, i.e. inside the brain
if ~isfield(sourcemodel, 'inside')
  if isfield(sourcemodel, 'tissue') && isfield(sourcemodel, 'tissuelabel')
    % this applies when basedoncentroids or movetocentroids
    % find the dipoles in the cortical or brain tissues
    cortex = find(ismember(headmodel.tissuelabel, {'gm', 'gray', 'brain'}));
    sourcemodel.inside = ismember(sourcemodel.tissue, cortex);
  else
    % this returns a boolean vector
    sourcemodel.inside = ft_inside_headmodel(sourcemodel.pos, headmodel, 'grad', sens, 'headshape', cfg.headshape, 'inwardshift', cfg.inwardshift);
  end
end % if inside

if strcmp(cfg.tight, 'yes')
  if ~isfield(sourcemodel, 'dim')
    ft_warning('cfg.tight only works for dipole positions on a regular 3D grid');
  else
    fprintf('%d dipoles inside, %d dipoles outside brain\n', sum(sourcemodel.inside), sum(~sourcemodel.inside));
    fprintf('making tight grid\n');
    boolvol = reshape(sourcemodel.inside, sourcemodel.dim);
    xsel    = squeeze(sum(sum(boolvol,3),2))>0;
    ysel    = squeeze(sum(sum(boolvol,3),1))>0;
    zsel    = squeeze(sum(sum(boolvol,2),1))>0;
    boolvol(xsel,ysel,zsel) = true; % update the volume to contain the to-be-selected entries
    sel     = boolvol(:);

    % update the boolean fields, this requires the original dim
    fn = booleanfields(sourcemodel);
    for i=1:numel(fn)
      sourcemodel.(fn{i}) = sourcemodel.(fn{i})(sel);
    end

    % update the grid locations that are marked as inside the brain
    sourcemodel.pos   = sourcemodel.pos(sel,:);
    sourcemodel.dim   = [sum(xsel) sum(ysel) sum(zsel)];
  end
end % if tight

fprintf('%d dipoles inside, %d dipoles outside brain\n', sum(sourcemodel.inside), sum(~sourcemodel.inside));

if istrue(cfg.checkinside) && ~any(sourcemodel.inside)
  ft_error('there are no dipoles inside the volume conductor')
end

% apply the symmetry constraint, i.e. add a symmetric dipole for each location that was defined sofar
if ~isempty(cfg.symmetry)
  if size(sourcemodel.pos,2)>3
    % sanity check, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3119
    ft_warning('the construction of a symmetric dipole model requires to start with a Nx3 description of the dipole positions, discarding subsequent columns');
    sourcemodel.pos = sourcemodel.pos(:,1:3);
  end
  if strcmp(cfg.symmetry, 'x')
    reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    mirror = [1 1 1 -1 1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 -x1 y
  elseif strcmp(cfg.symmetry, 'y')
    reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    mirror = [1 1 1 1 -1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 -y
  elseif strcmp(cfg.symmetry, 'z')
    reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    mirror = [1 1 1 1 1 -1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 y1
  else
    ft_error('unrecognized symmetry constraint');
  end
  fprintf('each source describes two dipoles with symmetry along %s axis\n', cfg.symmetry);
  % expand the number of parameters from one (3) to two dipoles (6)
  sourcemodel.pos = sourcemodel.pos(:,expand) .* repmat(mirror, size(sourcemodel.pos,1), 1);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble provenance sourcemodel
ft_postamble history    sourcemodel
ft_postamble savevar    sourcemodel


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function for basedonmri method to determine the inside
% returns a boolean vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inside = getinside(pos, mask)

% it might be that the box with the points does not completely fit into the mask
dim = size(mask);
sel = find(pos(:,1)<1 |  pos(:,1)>dim(1) | ...
  pos(:,2)<1 |  pos(:,2)>dim(2) | ...
  pos(:,3)<1 |  pos(:,3)>dim(3));
if isempty(sel)
  % use the efficient implementation
  inside = mask(sub2ind(dim, pos(:,1), pos(:,2), pos(:,3)));
else
  % only loop over the points that can be dealt with
  inside = false(size(pos,1), 1);
  for i=setdiff(1:size(pos,1), sel(:)')
    inside(i) = mask(pos(i,1), pos(i,2), pos(i,3));
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to return the fieldnames of the boolean fields in a
% segmentation, should work both for volumetric and for source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn = booleanfields(mri)

fn = fieldnames(mri);
isboolean = false(1,numel(fn));
for i=1:numel(fn)
  if islogical(mri.(fn{i})) && isequal(numel(mri.(fn{i})),prod(mri.dim))
    isboolean(i) = true;
  end
end
fn  = fn(isboolean);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute the centroids of the elements of a volumetric mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function centr = compute_centroids(headmodel)

% some of the fields can be copied over, fields that are specified but not present will be silently ignored
centr = keepfields(headmodel, {'tissue', 'tissuelabel', 'unit', 'coordsys'});

% the FEM model should have tetrahedrons or hexahedrons
if isfield(headmodel, 'tet')
  numtet = size(headmodel.tet, 1);
  fprintf('computing centroids for %d tetrahedrons\n', numtet);
  % compute the mean of the 4 corner points of the tetrahedrons
  centr.pos = squeeze(mean(reshape(headmodel.pos(headmodel.tet,:), numtet, 4, 3), 2));
elseif isfield(headmodel, 'hex')
  numhex = size(headmodel.hex, 1);
  fprintf('computing centroids for %d hexahedrons\n', numhex);
  % compute the mean of the 8 corner points of the hexahedrons
  centr.pos = squeeze(mean(reshape(headmodel.pos(headmodel.hex,:), numhex, 8, 3), 2));
else
  ft_error('the headmodel does not contain tetrahedrons or hexahedrons');
end
