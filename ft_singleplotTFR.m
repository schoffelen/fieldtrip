function [cfg] = ft_singleplotTFR(cfg, varargin)

% FT_SINGLEPLOTTFR plots the time-frequency representation of power of a
% single channel or the average over multiple channels.
%
% Use as
%   ft_singleplotTFR(cfg,data)
%
% The input freq structure should be a a time-frequency representation of
% power or coherence that was computed using the FT_FREQANALYSIS function.
%
% The configuration can have the following parameters:
%   cfg.parameter      = field to be plotted on z-axis, e.g. 'powspctrm' (default depends on data.dimord)
%   cfg.maskparameter  = field in the data to be used for masking of data, can be logical (e.g. significant data points) or numerical (e.g. t-values).
%                        (not possible for mean over multiple channels, or when input contains multiple subjects
%                        or trials)
%   cfg.maskstyle      = style used to masking, 'opacity', 'saturation', or 'outline' (default = 'opacity')
%                        'outline' can only be used with a logical cfg.maskparameter
%                        use 'saturation' or 'outline' when saving to vector-format (like *.eps) to avoid all sorts of image-problems
%   cfg.maskalpha      = alpha value between 0 (transparent) and 1 (opaque) used for masking areas dictated by cfg.maskparameter (default = 1)
%                        (will be ignored in case of numeric cfg.maskparameter or if cfg.maskstyle = 'outline')
%   cfg.masknans       = 'yes' or 'no' (default = 'yes')
%   cfg.xlim           = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.ylim           = 'maxmin' or [ymin ymax] (default = 'maxmin')
%   cfg.zlim           = plotting limits for color dimension, 'maxmin', 'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
%   cfg.baseline       = 'yes', 'no' or [time1 time2] (default = 'no'), see FT_FREQBASELINE
%   cfg.baselinetype   = 'absolute', 'relative', 'relchange', 'normchange', 'db' or 'zscore' (default = 'absolute')
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                        see FT_CHANNELSELECTION for details
%   cfg.title          = string, title of plot
%   cfg.refchannel     = name of reference channel for visualising connectivity, can be 'gui'
%   cfg.fontsize       = font size of title (default = 8)
%   cfg.hotkeys        = enables hotkeys (leftarrow/rightarrow/uparrow/downarrow/pageup/pagedown/m) for dynamic zoom and translation (ctrl+) of the axes and color limits
%   cfg.colormap       = string, or Nx3 matrix, see FT_COLORMAP
%   cfg.colorbar       = 'yes', 'no' (default = 'yes')
%   cfg.colorbartext   = string indicating the text next to colorbar
%   cfg.interactive    = interactive plot 'yes' or 'no' (default = 'yes')
%                        In a interactive plot you can select areas and produce a new
%                        interactive plot when a selected area is clicked. Multiple areas
%                        can be selected by holding down the SHIFT key.
%   cfg.position       = location and size of the figure, specified as [left bottom width height] (default is automatic)
%   cfg.renderer       = string, 'opengl', 'zbuffer', 'painters', see RENDERERINFO (default is automatic, try 'painters' when it crashes)
%   cfg.directionality = '', 'inflow' or 'outflow' specifies for
%                       connectivity measures whether the inflow into a
%                       node, or the outflow from a node is plotted. The
%                       (default) behavior of this option depends on the dimor
%                       of the input data (see below).
%   cfg.figure         = 'yes', 'no', or 'subplot',  whether to open a new figure. You can also specify a figure
%                        handle from FIGURE, GCF or SUBPLOT. (default = 'yes'). With multiple data inputs, 'subplot'
%                        will make subplots in a single figure.
%
% The following options for the scaling of the EEG, EOG, ECG, EMG, MEG and NIRS channels
% is optional and can be used to bring the absolute numbers of the different
% channel types in the same range (e.g. fT and uV). The channel types are determined
% from the input data using FT_CHANNELSELECTION.
%   cfg.eegscale       = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale       = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale       = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale       = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale       = number, scaling to apply to the MEG channels prior to display
%   cfg.gradscale      = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.magscale       = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.nirsscale      = number, scaling to apply to the NIRS channels prior to display
%   cfg.mychanscale    = number, scaling to apply to the channels specified in cfg.mychan
%   cfg.mychan         = Nx1 cell-array with selection of channels
%   cfg.chanscale      = Nx1 vector with scaling factors, one per channel specified in cfg.channel
%
% For the plotting of directional connectivity data the cfg.directionality option determines what is plotted. The default
% value and the supported functionality depend on the dimord of the input data. If the input data is of dimord 'chan_chan_XXX',
% the value of directionality determines whether, given the reference channel(s), the columns (inflow), or rows (outflow) are
% selected for plotting. In this situation the default is 'inflow'. Note that for undirected measures, inflow and outflow should
% give the same output. If the input data is of dimord 'chancmb_XXX', the value of directionality determines whether the rows in
% data.labelcmb are selected. With 'inflow' the rows are selected if the refchannel(s) occur in the right column, with 'outflow'
% the rows are selected if the refchannel(s) occur in the left column of the labelcmb-field. Default in this case is '', which
% means that all rows are selected in which the refchannel(s) occur. This is to robustly support linearly indexed undirected
% connectivity metrics. In the situation where undirected connectivity measures are linearly indexed, specifying 'inflow' or 
% outflow' can result in unexpected behavior.
%
% See also FT_SINGLEPLOTER, FT_MULTIPLOTER, FT_MULTIPLOTTFR, FT_TOPOPLOTER, FT_TOPOPLOTTFR

% Copyright (C) 2005-2025, F.C. Donders Centre
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPERS NOTE: This code is organized in a similar fashion for multiplot/singleplot/topoplot
% and for ER/TFR and should remain consistent over those 6 functions.
% Section 1: general cfg handling that is independent from the data
% Section 2: data handling, this also includes converting bivariate (chan_chan and chancmb) into univariate data
% Section 3: select the data to be plotted and determine min/max range
% Section 4: do the actual plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section 1: general cfg handling that is independent from the data

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    varargin
ft_preamble provenance varargin

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
Ndata = numel(varargin);
for i=1:Ndata
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'freq');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',   {'channels', 'trial'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'unused',      {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'matrixside',     'directionality'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'zlim', 'absmax', 'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'directionality', 'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'directionality', 'feedback',    'inflow'});
cfg = ft_checkconfig(cfg, 'renamed',     {'channelindex',   'channel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'channelname',    'channel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'cohrefchannel',  'refchannel'});
cfg = ft_checkconfig(cfg, 'renamed',	   {'zparam',         'parameter'});
cfg = ft_checkconfig(cfg, 'renamed',     {'newfigure',      'figure'});

% Set the defaults
cfg.baseline       = ft_getopt(cfg, 'baseline',      'no');
cfg.baselinetype   = ft_getopt(cfg, 'baselinetype',  'absolute');
cfg.trials         = ft_getopt(cfg, 'trials',        'all', 1);
cfg.xlim           = ft_getopt(cfg, 'xlim',          'maxmin');
cfg.ylim           = ft_getopt(cfg, 'ylim',          'maxmin');
cfg.zlim           = ft_getopt(cfg, 'zlim',          'maxmin');
cfg.fontsize       = ft_getopt(cfg, 'fontsize',       8);
cfg.interpreter    = ft_getopt(cfg, 'interpreter',   'none');
cfg.colorbar       = ft_getopt(cfg, 'colorbar',      'yes');
cfg.colormap       = ft_getopt(cfg, 'colormap',       'default');
cfg.colorbartext   = ft_getopt(cfg, 'colorbartext',  '');
cfg.interactive    = ft_getopt(cfg, 'interactive',   'yes');
cfg.hotkeys        = ft_getopt(cfg, 'hotkeys',       'yes');
cfg.maskalpha      = ft_getopt(cfg, 'maskalpha',      1);
cfg.maskparameter  = ft_getopt(cfg, 'maskparameter',  []);
cfg.maskstyle      = ft_getopt(cfg, 'maskstyle',     'opacity');
cfg.channel        = ft_getopt(cfg, 'channel',       'all');
cfg.title          = ft_getopt(cfg, 'title',          []);
cfg.masknans       = ft_getopt(cfg, 'masknans',      'yes');
cfg.directionality = ft_getopt(cfg, 'directionality', []);
cfg.figurename     = ft_getopt(cfg, 'figurename',     []);
cfg.parameter      = ft_getopt(cfg, 'parameter',     'powspctrm');
cfg.renderer       = ft_getopt(cfg, 'renderer',       []); % let MATLAB decide on the default
cfg.figure         = ft_getopt(cfg, 'figure',         'yes');

% this is needed for the figure title
if isfield(cfg, 'dataname') && ~isempty(cfg.dataname)
  dataname = cfg.dataname;
elseif isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  dataname = cfg.inputfile;
elseif nargin>1
  dataname = arrayfun(@inputname, 2:nargin, 'UniformOutput', false);
else
  dataname = {};
end

makesubplots = false;
if Ndata==1 && isequal(cfg.figure, 'subplot')
  % overrule this setting
  cfg.figure = 'yes';
elseif Ndata>1 && isequal(cfg.figure, 'subplot')
  makesubplots = true;
end
  
%% Section 2: data handling, this also includes converting bivariate (chan_chan and chancmb) into univariate data

hastime = isfield(varargin{1}, 'time');
hasfreq = isfield(varargin{1}, 'freq');

assert((hastime && hasfreq), 'please use ft_singleplotER for time-only or frequency-only data');

xparam = ft_getopt(cfg, 'xparam', 'time');
yparam = ft_getopt(cfg, 'yparam', 'freq');

% check whether rpt/subj is present and remove if necessary
dimord = getdimord(varargin{1}, cfg.parameter);
dimtok = tokenize(dimord, '_');
hasrpt = any(ismember(dimtok, {'rpt' 'subj'}));

if ~hasrpt
  assert(isequal(cfg.trials, 'all') || isequal(cfg.trials, 1), 'incorrect specification of cfg.trials for data without repetitions');
else
  assert(~isempty(cfg.trials), 'empty specification of cfg.trials for data with repetitions');
end

% parse cfg.channel
if isfield(cfg, 'channel') && isfield(varargin{1}, 'label')
  cfg.channel = ft_channelselection(cfg.channel, varargin{1}.label);
elseif isfield(cfg, 'channel') && isfield(varargin{1}, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(varargin{1}.labelcmb(:)));
end

% Apply baseline correction:
if ~strcmp(cfg.baseline, 'no')
  for i=1:Ndata
    tmpcfg = keepfields(cfg, {'baseline', 'baselinetype', 'baselinewindow', 'demean', 'parameter', 'channel'});
    % keep mask-parameter if it is set
    if ~isempty(cfg.maskparameter)
      tempmask = varargin{i}.(cfg.maskparameter);
    end
    varargin{i} = ft_freqbaseline(tmpcfg, varargin{i});
    % put mask-parameter back if it is set
    if ~isempty(cfg.maskparameter)
      varargin{i}.(cfg.maskparameter) = tempmask;
    end
  end
end

% channels should NOT be selected and averaged here, since a topoplot might follow in interactive mode
tmpcfg = keepfields(cfg, {'trials', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
if hasrpt
  tmpcfg.avgoverrpt = 'yes';
else
  tmpcfg.avgoverrpt = 'no';
end
tmpvar = varargin{1};
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information, don't keep the ft_selectdata details
[tmpcfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

if isfield(tmpvar, cfg.maskparameter) && ~isfield(varargin{1}, cfg.maskparameter)
  % the mask parameter is not present after ft_selectdata, because it is
  % not included in all input arguments. Make the same selection and copy
  % it over
  tmpvar = ft_selectdata(tmpcfg, tmpvar);
  varargin{1}.(cfg.maskparameter) = tmpvar.(cfg.maskparameter);
end

clear tmpvar tmpcfg dimord dimtok hastime hasfreq hasrpt

% ensure that the preproc specific options are located in the cfg.preproc
% substructure, but also ensure that the field 'refchannel' remains at the
% highest level in the structure. This is a little hack by JM because the field
% refchannel can relate to connectivity or to an EEg reference.

if isfield(cfg, 'refchannel'), refchannelincfg = cfg.refchannel; cfg = rmfield(cfg, 'refchannel'); end
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});
if exist('refchannelincfg', 'var'), cfg.refchannel  = refchannelincfg; end

if ~isempty(cfg.preproc)
  % preprocess the data, i.e. apply filtering, baselinecorrection, etc.
  fprintf('applying preprocessing options\n');
  if ~isfield(cfg.preproc, 'feedback')
    cfg.preproc.feedback = cfg.interactive;
  end
  for i=1:Ndata
    varargin{i} = ft_preprocessing(cfg.preproc, varargin{i});
  end
end

% Handle the bivariate case
dimord = getdimord(varargin{1}, cfg.parameter);
if startsWith(dimord, 'chan_chan_') || startsWith(dimord, 'chancmb_')
  % convert the bivariate data to univariate and call this plotting function again
  cfg.originalfunction = 'ft_singleplotTFR';
  cfg.trials = 'all'; % trial selection has been taken care off
  bivariate_common(cfg, varargin{:});
  return
end

% Apply channel-type specific scaling
fn = fieldnames(cfg);
fn = setdiff(fn, {'skipscale', 'showscale', 'gridscale'}); % these are for the layout and plotting, not for CHANSCALE_COMMON
fn = fn(endsWith(fn, 'scale') | startsWith(fn, 'mychan') | strcmp(fn, 'channel') | strcmp(fn, 'parameter'));
tmpcfg = keepfields(cfg, fn);
if ~isempty(tmpcfg)
  for i=1:Ndata
    varargin{i} = chanscale_common(tmpcfg, varargin{i});
  end
  % remove the scaling fields from the configuration, to prevent them from being called again in interactive mode
  % but keep the parameter and channel field
  cfg = removefields(cfg, setdiff(fn, {'parameter', 'channel'}));
else
  % do nothing
end

%% Section 3: select the data to be plotted and determine min/max range

% Take the desided subselection of channels, this is the same in all datasets
[selchan] = match_str(varargin{1}.label, cfg.channel);

% Get physical min/max range of x, i.e. time
if strcmp(cfg.xlim, 'maxmin')
  % Find maxmin throughout all varargins:
  xmin = [];
  xmax = [];
  for i=1:Ndata
    xmin = min([xmin varargin{i}.(xparam)]);
    xmax = max([xmax varargin{i}.(xparam)]);
  end
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Get the index of the nearest bin, this is the same in all datasets
xminindx = nearest(varargin{1}.(xparam), xmin);
xmaxindx = nearest(varargin{1}.(xparam), xmax);
xmin = varargin{1}.(xparam)(xminindx);
xmax = varargin{1}.(xparam)(xmaxindx);
selx = xminindx:xmaxindx;
xval = varargin{1}.(xparam)(selx);

% Get physical min/max range of y, i.e. frequency
if strcmp(cfg.ylim, 'maxmin')
  % Find maxmin throughout all varargins:
  ymin = [];
  ymax = [];
  for i=1:Ndata
    ymin = min([ymin varargin{i}.(yparam)]);
    ymax = max([ymax varargin{i}.(yparam)]);
  end
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% Get the index of the nearest bin
yminindx = nearest(varargin{1}.(yparam), ymin);
ymaxindx = nearest(varargin{1}.(yparam), ymax);
ymin = varargin{1}.(yparam)(yminindx);
ymax = varargin{1}.(yparam)(ymaxindx);
sely = yminindx:ymaxindx;
yval = varargin{1}.(yparam)(sely);

% test if X and Y are linearly spaced (to within 10^-12): % FROM UIMAGE
dx = min(diff(xval));  % smallest interval for X
dy = min(diff(yval));  % smallest interval for Y
evenx = all(abs(diff(xval)/dx-1)<1e-12);     % true if X is linearly spaced
eveny = all(abs(diff(yval)/dy-1)<1e-12);     % true if Y is linearly spaced

if ~evenx || ~eveny
  ft_warning('(one of the) axis is/are not evenly spaced, but plots are made as if axis are linear')
end

% masking is only possible for evenly spaced axis
if strcmp(cfg.masknans, 'yes') && (~evenx || ~eveny)
  ft_warning('(one of the) axis are not evenly spaced -> nans cannot be masked out -> cfg.masknans is set to ''no'';')
  cfg.masknans = 'no';
end

% the usual data is chan_freq_time, but other dimords should also work
dimtok = tokenize(dimord, '_');
for i=1:Ndata
  data = varargin{i};

  datamatrix = data.(cfg.parameter);
  [c, ia, ib] = intersect({'chan', yparam, xparam}, dimtok, 'stable');
  datamatrix = permute(datamatrix, ib);
  datamatrix = datamatrix(selchan, sely, selx);

  if ~isempty(cfg.maskparameter) && isfield(data, cfg.maskparameter)
    maskmatrix = data.(cfg.maskparameter)(selchan, sely, selx);
    if islogical(maskmatrix) && any(strcmp(cfg.maskstyle, {'saturation', 'opacity'}))
      maskmatrix = double(maskmatrix);
      maskmatrix(~maskmatrix) = cfg.maskalpha;
    elseif isnumeric(maskmatrix)
      if strcmp(cfg.maskstyle, 'outline')
        ft_error('Outline masking with a numeric cfg.maskparameter is not supported. Please use a logical mask instead.')
      end
      if cfg.maskalpha ~= 1
        ft_warning('Using field "%s" for masking, cfg.maskalpha is ignored.', cfg.maskparameter)
      end
      % scale mask between 0 and 1
      minval = min(maskmatrix(:));
      maxval = max(maskmatrix(:));
      maskmatrix = (maskmatrix - minval) / (maxval-minval);
    end
  else
    % create an Nx0x0 matrix
    maskmatrix = zeros(length(selchan), 0, 0);
  end

  %% Section 4: do the actual plotting
  if makesubplots
    % make multiple plots in a single figure
    nyplot = ceil(sqrt(Ndata));
    nxplot = ceil(Ndata./nyplot);
    cfg.figure = subplot(nxplot, nyplot, i);
  end
  
  % open a new figure, or add it to the existing one
  % note that in general adding a TFR to an existing one does not make sense, since they will overlap
  open_figure(keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'}));

  zval = mean(datamatrix, 1); % over channels
  zval = reshape(zval, size(zval,2), size(zval,3));
  mask = squeeze(mean(maskmatrix, 1)); % over channels

  % Get physical z-axis range (color axis):
  if strcmp(cfg.zlim, 'maxmin')
    zmin = min(zval(:), [], 'omitnan');
    zmax = max(zval(:), [], 'omitnan');
  elseif strcmp(cfg.zlim, 'maxabs')
    zmin = -max(abs(zval(:)), [], 'omitnan');
    zmax =  max(abs(zval(:)), [], 'omitnan');
  elseif strcmp(cfg.zlim, 'zeromax')
    zmin = 0;
    zmax = max(zval(:), [], 'omitnan');
  elseif strcmp(cfg.zlim, 'minzero')
    zmin = min(zval(:), [], 'omitnan');
    zmax = 0;
  else
    zmin = cfg.zlim(1);
    zmax = cfg.zlim(2);
  end

  % Draw the data and mask NaN's if requested
  if isequal(cfg.masknans, 'yes') && isempty(cfg.maskparameter)
    nans_mask = ~isnan(zval);
    mask = double(nans_mask);
    ft_plot_matrix(xval, yval, zval, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask)
  elseif isequal(cfg.masknans, 'yes') && ~isempty(cfg.maskparameter)
    nans_mask = ~isnan(zval);
    mask = mask .* nans_mask;
    mask = double(mask);
    ft_plot_matrix(xval, yval, zval, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask)
  elseif isequal(cfg.masknans, 'no') && ~isempty(cfg.maskparameter)
    mask = double(mask);
    ft_plot_matrix(xval, yval, zval, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask)
  else
    ft_plot_matrix(xval, yval, zval, 'clim', [zmin zmax], 'tag', 'cip')
  end

  % check if the colormap is in the proper format and set it
  if ~isequal(cfg.colormap, 'default')
    if ischar(cfg.colormap)
      cfg.colormap = ft_colormap(cfg.colormap);
    elseif iscell(cfg.colormap)
      cfg.colormap = ft_colormap(cfg.colormap{:});
    elseif isnumeric(cfg.colormap) && size(cfg.colormap,2)~=3
      ft_error('colormap must be a Nx3 matrix');
    end
    set(gcf, 'colormap', cfg.colormap);
  end

  axis xy

  if isequal(cfg.colorbar, 'yes')
    c = colorbar;
    ylabel(c, cfg.colorbartext);
  end

  % Set callback to adjust color axis
  if strcmp('yes', cfg.hotkeys)
    %  Attach data and cfg to figure and attach a key listener to the figure
    set(gcf, 'KeyPressFcn', {@key_sub, xmin, xmax, ymin, ymax, zmin, zmax})
  end

  % Create axis title containing channel name(s) and channel number(s):
  if ~isempty(cfg.title)
    t = cfg.title;
  else
    if length(cfg.channel) == 1
      t = [char(cfg.channel) ' / ' num2str(selchan) ];
    else
      t = sprintf('mean(%0s)', join_str(', ', cfg.channel));
    end
  end
  title(t, 'fontsize', cfg.fontsize, 'interpreter', cfg.interpreter);

  % set the figure window title, add channel labels if number is small
  if isempty(get(gcf, 'Name'))
    if length(selchan) < 5
      chans = join_str(', ', cfg.channel);
    else
      chans = '<multiple channels>';
    end
    if ~isempty(cfg.figurename)
      set(gcf, 'name', cfg.figurename);
      set(gcf, 'NumberTitle', 'off');
    elseif ~isempty(dataname)
      set(gcf, 'Name', sprintf('%d: %s: %s (%s)', double(gcf), mfilename, join_str(', ', dataname), chans));
      set(gcf, 'NumberTitle', 'off');
    else
      set(gcf, 'Name', sprintf('%d: %s (%s)', double(gcf), mfilename, chans));
      set(gcf, 'NumberTitle', 'off');
    end
  end

  axis tight

  % Make the figure interactive
  if strcmp(cfg.interactive, 'yes')
    % add the cfg/data information to the figure under identifier linked to this axis
    ident             = ['axh' num2str(round(sum(clock.*1e6)))]; % unique identifier for this axis
    set(gca, 'tag',ident);

    % ensure that the function that is called knows about the subplot setting
    if makesubplots
      cfg.figure = 'subplot';
    end
    info                  = guidata(gcf);
    info.(ident).dataname = dataname;
    info.(ident).cfg      = cfg;
    info.(ident).varargin = varargin;
    guidata(gcf, info);
    set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR}, 'event', 'WindowButtonMotionFcn'});
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous data
ft_postamble provenance
ft_postamble savefig

% add a menu to the figure, but only if the current figure does not have subplots
menu_fieldtrip(gcf, cfg, false);

if ~ft_nargout
  % don't return anything
  clear cfg
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting a time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotTFR(range, varargin)
% fetch cfg/data based on axis indentifier given as tag
ident  = get(gca, 'tag');
info   = guidata(gcf);
cfg    = info.(ident).cfg;
varargin = info.(ident).varargin;
if ~isempty(range)
  cfg = removefields(cfg, 'inputfile');   % the reading has already been done and varargin contains the data
  cfg = removefields(cfg, 'showlabels');  % this is not allowed in topoplotER
  cfg.trials = 'all';                     % trial selection has already been taken care of
  cfg.baseline = 'no';                    % make sure the next function does not apply a baseline correction again
  cfg.channel = 'all';                    % make sure the topo displays all channels, not just the ones in this singleplot
  cfg.comment = 'auto';
  cfg.dataname = info.(ident).dataname;   % put data name in here, this cannot be resolved by other means
  cfg.xlim = range(1:2);
  cfg.ylim = range(3:4);
  fprintf('selected cfg.xlim = [%f %f]\n', cfg.xlim(1), cfg.xlim(2));
  fprintf('selected cfg.ylim = [%f %f]\n', cfg.ylim(1), cfg.ylim(2));
  % ensure that the new figure appears at the same position
  cfg.position = get(gcf, 'Position');
  if isfield(cfg, 'figure') && isequal(cfg.figure, 'subplot')
    figure('position', cfg.position);
  else
    cfg.figure = 'yes';
  end
  ft_topoplotTFR(cfg, varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata, varargin)
xlimits = xlim;
ylimits = ylim;
climits = clim;
incr_x = abs(xlimits(2) - xlimits(1)) /10;
incr_y = abs(ylimits(2) - ylimits(1)) /10;
incr_c = abs(climits(2) - climits(1)) /10;

if length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:}, 'control')
  % TRANSLATE by 10%
  switch eventdata.Key
    case 'pageup'
      clim([min(clim)+incr_c max(clim)+incr_c]);
    case 'pagedown'
      clim([min(clim)-incr_c max(clim)-incr_c]);
    case 'leftarrow'
      xlim([xlimits(1)+incr_x xlimits(2)+incr_x])
    case 'rightarrow'
      xlim([xlimits(1)-incr_x xlimits(2)-incr_x])
    case 'uparrow'
      ylim([ylimits(1)-incr_y ylimits(2)-incr_y])
    case 'downarrow'
      ylim([ylimits(1)+incr_y ylimits(2)+incr_y])
  end % switch
else
  % ZOOM by 10%
  switch eventdata.Key
    case 'pageup'
      clim([min(clim)-incr_c max(clim)+incr_c]);
    case 'pagedown'
      clim([min(clim)+incr_c max(clim)-incr_c]);
    case 'leftarrow'
      xlim([xlimits(1)-incr_x xlimits(2)+incr_x])
    case 'rightarrow'
      xlim([xlimits(1)+incr_x xlimits(2)-incr_x])
    case 'uparrow'
      ylim([ylimits(1)-incr_y ylimits(2)+incr_y])
    case 'downarrow'
      ylim([ylimits(1)+incr_y ylimits(2)-incr_y])
    case 'm'
      xlim([varargin{1} varargin{2}])
      ylim([varargin{3} varargin{4}])
      clim([varargin{5} varargin{6}]);
  end % switch
end % if
