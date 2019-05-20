function [output] = ft_connectivity_mim(input, varargin)

% FT_CONNECTIVITY_MIM computes the multivariate interaction measure from a
% data-matrix containing the cross-spectral density. This implements the method
% described in Ewald et al., Estimating true brain connectivity from EEG/MEG data
% invariant to linear and static trasformations in sensor space. Neuroimage, 2012;
% 476:488.
%
% Use as
%   [m] = hcp_connectivity_mim(input, ...)
%
% The input data should be an array organized as
%   Channel x Channel x Frequency
%
% Additional optional input arguments come as key-value pairs:
%   indices   = 1xN vector with indices of the groups to which the channels belong,
%               e.g. [1 1 2 2] for a 2-by-2 connectivity between planar MEG channels.
%
% The output m contains the Channel*Channel connectivity measure.
%
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2011-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
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

indices = ft_getopt(varargin, 'indices');

if isempty(indices) && isequal(size(input), [2 2])
  % simply assume two channels
  indices = [1 2];
else
  
end

siz = [size(input) 1];
if siz(1)~=siz(2)
  error('currently the input is expected to be of size NxNxM');
end

uindices = unique(indices(:));
N = numel(uindices);
if ~isequal(uindices,(1:N)')
  error('indices values should range from 1 to the number of channels/voxels');
end

input_r = real(input);
input_i = imag(input);

output = zeros([N N siz(3:end)]);
for k = 1:prod(siz(3:end))
  for m1 = 1:N
    for m2 = 1:N
      indx1 = indices==m1;
      indx2 = indices==m2;
      cs_aa_re = input_r(indx1,indx1,k);
      cs_bb_re = input_r(indx2,indx2,k);
      cs_ab_im = input_i(indx1,indx2,k);
      output(m1,m2,k) = trace((cs_aa_re\cs_ab_im/cs_bb_re)*transpose(cs_ab_im));
    end
  end
end
