function [p] = ft_connectivity_plm(input, varargin)

% FT_CONNECTIVITY_PLM computes the phase linearity measurement from a cell
% array of time-domain data, where each cell is an epoch
%
% Use as
%   [p] = ft_connectivity_plm(input, ...)
%
% The input data input should be organized as a cell-array of nchan x ntime signals
%
% Additional optional input arguments come as key-value pairs:
%   bandwidth	=	scalar, half-bandwidth parameter: the frequency range
%			across which to integrate
%   fsample     =       sampling frequency, needed to convert bandwidth to number of bins
%
% The output p contains the phase linearity measurement in the [0, 1] interval.
% The output p is organized as a 3D matrix of nchan x  nchan x ntime dimensions.
%
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2018, Fabio Baselice, Pierpaolo Sorrentino, Jan-Mathijs Schoffelen
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

% the sequence of steps is as follows:
%  - Hilbert transformation
%  - multiply with complex conjugate
%  - fft
%  - convert bandwidth parameter to number of bins
%  - integrate over bandwidth


% NOTE BY JM: if the user inputs data with different length trials, the fft per trial is going 
% to have different frequency resolutions, which is not good. Better to throw an error in that 
% case.
nsmp = cellfun('size', input, 2);
assert(all(nsmp==nsmp(1)), 'currently there is no support for input, where the trials are of different length'); 

ntime=numel(input);
for k = 1:numel(input)
  input{k} = hilbert(input{k}')';
end

nchan=size(input{1},1);
trial_length=size(input{1},2);
ph_min=0.1;        % Eps of Eq.(17) of the manuscript
f=(fs/trial_length)*(0:(trial_length-1));
f_integr=(abs(f)<B) | (abs(f-fs)<B);
p=zeros(nchan, nchan, ntime);

for ktime=1:ntime
    for kchan1=1:(nchan-1)
        for kchan2=(kchan1+1):nchan
            temp=fft(input{ktime}(kchan1,:).*conj(input{ktime}(kchan2,:)));    % NOTE BY FB: The inner cycle can be vectorized
            temp(1)=temp(1).*(abs(angle(temp(1)))>ph_min);  % Volume conduction suppression
            temp=(abs(temp)).^2;
            p_temp=sum(temp(f_integr))./sum(temp);
            p(kchan1, kchan2, ktime)=p_temp;
            p(kchan2, kchan1, ktime)=p_temp;
        end
    end
end

