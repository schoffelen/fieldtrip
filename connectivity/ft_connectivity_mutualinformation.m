function output = ft_connectivity_mutualinformation(input, varargin)

% FT_CONNECTIVITY_MUTUALINFORMATION computes mutual information using
% either the information breakdown toolbox (ibtb), as described in Magri
% et al., BMC Neuroscience 2009, 1471-2202, or Robin Ince's Gaussian copula
% based parametric approach (gcmi). The function is a helper function for
% FT_CONNECTIVITYANALYSIS. As a standalone function, it could be used as
% follows:
%
% mi = ft_connectivity_mutualinformation(data, varargin)
%
% The input data is a Nchan x Nobservations matrix.
%
% Additional input arguments come as key-value pairs:
%   method     = string, 'ibtb' (default), or 'gcmi'.
%
% The following arguments pertain to the 'ibtb' method:
%   histmethod = The way that histograms are generated from the data. Possible values
%                are 'eqpop' (default), 'eqspace', 'ceqspace', 'gseqspace'.
%                See the help of the 'binr' function in the ibtb toolbox for more information.
%   numbin     = scalar value. The number of bins used to create the histograms needed for
%                the entropy computations
%   opts       = structure that is passed on to the 'information' function in the ibtb
%                toolbox. See the help of that function for more information.
%   refindx    = scalar value or 'all'. The channel that is used as 'reference channel'.
%
% The output contains the estimated mutual information between all channels and the reference channel(s).

% Copyright (C) 2016 Donders Institute, Jan-Mathijs Schoffelen
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

method  = ft_getopt(varargin, 'method',  'ibtb'); % can be gcmi
refindx = ft_getopt(varargin, 'refindx', 'all', 1);
lags    = ft_getopt(varargin, 'lags',    0);  % shift of data w.r.t. reference, in samples
tra     = ft_getopt(varargin, 'tra',     []); % 1/0-matrix for multivariate combination Nnew x Norg, where Norg = size(input,1)

% check whether the combined options work out
if ~isempty(tra)
  tra = full(tra)>0;
  if strcmp(method, 'ibtb') && ~isequal(tra,eye(size(tra,1))>0),
    error('method ''ibtb'' in combination with a non-identity ''tra'' is not possible');
  end
else
  tra = eye(size(input,1))>0;
end

% ensure that the refindx is numeric, defaults to 1:size(input,1), i.e. do
% all-to-all
if (ischar(refindx) && strcmp(refindx, 'all')) || isempty(refindx)
  refindx = (1:size(tra,1))';
end

% do not allow anything else than a scalar, or 1:nchan as refindx
if numel(refindx)~=1 && numel(refindx)~=size(tra,1)
  error('mi can only be computed using a single, or all channels as reference');
end

switch method
  case 'ibtb'
    % check whether the required toolbox is available
    ft_hastoolbox('ibtb', 1);
    
    % set some options
    histmethod = ft_getopt(varargin, 'histmethod', 'eqpop');
    numbin     = ft_getopt(varargin, 'numbin',     10);
    
    % set some additional options that pertain to the algorithmic details of the
    % mutual information computation, see the documentation of ibtb
    opts        = ft_getopt(varargin, 'opts', []);
    opts.nt     = ft_getopt(opts, 'nt', []);
    opts.method = ft_getopt(opts, 'method', 'dr');
    opts.bias   = ft_getopt(opts, 'bias',   'pt');
    
    % deal with NaNs in the input data, e.g. trial boundaries
    finitevals = isfinite(input);
    
    nchans = size(tra,1); 
    n      = size(input, 2);
    output = zeros(nchans, numel(refindx), numel(lags)) + nan;
    
    % for each lag
    for m = 1:numel(lags)
      fprintf('computing mutualinformation for time lag in samples %d\n', lags(m));
      
      % get the samples for the relative shifts
      beg1 = max(0, lags(m))  + 1;
      beg2 = max(0, -lags(m)) + 1;
      n1   = n-abs(lags(m));
        
      end1 = beg1+n1-1;
      end2 = beg2+n1-1;
      
      for p = 1:numel(refindx)
        tmprefdata = nan(sum(tra(refindx(p),:)),n);
        tmprefdata(:, beg1:end1) = input(tra(refindx(p),:), beg2:end2);
        
        finitevals2 = sum(finitevals,1)&sum(isfinite(tmprefdata),1); % this conservatively takes only the non-nan samples across all input data channels
        
        tmpinput    = input(:,finitevals2);
        tmprefdata  = tmprefdata(:,finitevals2);
      
        % discretize signal1
        tmprefdata = binr(tmprefdata, sum(finitevals2), numbin, histmethod);
      
        for k = setdiff(1:size(tmpinput,1),refindx(p))
          signal2 = tmpinput(k,:);
        
          % represent signal2 in bins according to signal1's discretization
          R = zeros(1,3,numbin);
          for j = 1:numbin
            nr         = tmprefdata==j-1;
            opts.nt(j) = sum(nr);
            R(1, 1:opts.nt(j),j) = signal2(nr);
          end
        
          % discretize signal2 and compute mi
          R2 = binr(R, opts.nt', numbin, histmethod);
          output(k,p,m) = information(R2, opts, 'I'); % this computes mutual information
        end
      end
    end
    
  case 'gcmi'
    ft_hastoolbox('gcmi', 1);
    
    % set some options
    cmplx = ft_getopt(varargin, 'complex', 'complex'); % this is only used if data are complex-valued

    % deal with NaNs in the input data, e.g. trial boundaries
    finitevals = isfinite(input);
    
    % verify whether data is complex-valued, check the inputs, and adjust
    % the input data
    if ~all(imag(input(:))==0),
      % a tra deviating from I is currently not supported: ask Robin how to
      % deal with this, if possible at all
      if ~isequal(tra,eye(size(tra,1))>0),
        error('complex-valued input data in combination with multivariate signals is not supported');
      end
      switch cmplx
        case 'complex'
          % tease apart the real/imag parts, treat as 2D-variable, and
          % ensure the nans to behave
          input(~finitevals) = nan+1i.*nan;
          input = cat(1, real(input), imag(input));
          tra   = cat(2, tra, tra);
          finitevals = cat(1, finitevals, finitevals);
        case 'abs'
          % take the amplitude
          input = abs(input);
        case 'angle'
          % tease apart the real/imag parts, after amplitude normalization,
          % and ensure the nans to behave
          input(~finitevals) = nan+1i.*nan;
          input = input./abs(input);
          input = cat(1, real(input), imag(input));
          tra   = cat(2, tra, tra);
          finitevals = cat(1, finitevals, finitevals);
        otherwise
          error('unsupported value for ''complex''');
      end
    end
  
    nchans = size(tra,1); 
    n      = size(input, 2);
    output = zeros(nchans, numel(refindx), numel(lags)) + nan;
    
    % for each lag
    for m = 1:numel(lags)
      fprintf('computing mutualinformation for time lag in samples %d\n', lags(m));
      
      % get the samples for the relative shifts
      beg1 = max(0, lags(m))  + 1;
      beg2 = max(0, -lags(m)) + 1;
      n1   = n-abs(lags(m));
        
      end1 = beg1+n1-1;
      end2 = beg2+n1-1;
      
      for p = 1:numel(refindx)
        tmprefdata = nan(sum(tra(refindx(p),:)),n);
        tmprefdata(:, beg1:end1) = input(tra(refindx(p),:), beg2:end2);
        
        finitevals2 = sum(finitevals,1)&sum(isfinite(tmprefdata),1); % this conservatively takes only the non-nan samples across all input data channels
        
        tmpinput    = copnorm(input(:,finitevals2)')';
        tmprefdata  = copnorm(tmprefdata(:,finitevals2)')';
        for k = setdiff(1:size(tra,1),refindx(p))
          output(k,p,m) = mi_gg(tmpinput(tra(k,:),:)',tmprefdata');
        end
      end
    end
  otherwise
end

if numel(refindx)==1,
  siz    = [size(output) 1];
  output = reshape(output,[siz([1 3])]);
end
