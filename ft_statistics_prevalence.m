function [stat, cfg] = ft_statistics_prevalence(cfg, dat, design, varargin)

% FT_STATISTICS_PREVALENCE performs nonparametric prevalence inference
% based on:
% C Allefeld, K Görgen and JD Haynes
% Valid population inference for information-based imaging: From the
% second-level t-test to prevalence inference
%
% Monte-Carlo estimates of the significance probabilities and/or critical values
% from the permutation distribution. This function should not be called
% directly, instead you should call the function that is associated with the
% type of data on which you want to perform the test.
%
% Use as
%   stat = ft_timelockstatistics(cfg, data1, data2, data3, ...)
%   stat = ft_freqstatistics    (cfg, data1, data2, data3, ...)
%   stat = ft_sourcestatistics  (cfg, data1, data2, data3, ...)
%
% Where the data is obtained from FT_TIMELOCKANALYSIS, FT_FREQANALYSIS
% or FT_SOURCEANALYSIS respectively, or from FT_TIMELOCKGRANDAVERAGE,
% FT_FREQGRANDAVERAGE or FT_SOURCEGRANDAVERAGE respectively and with
% cfg.method = 'prevalence'
%
% The configuration options that can be specified are:
%   cfg.numrandomization1 = number of 1st level permutations
%   cfg.numrandomization2 = number of 2nd level permutations
%   cfg.alpha            = number, critical value for rejecting the null-hypothesis per tail (default = 0.05)
%   cfg.ivar             = number or list with indices, independent variable(s)
%   cfg.uvar             = number or list with indices, unit variable(s) (subjects)
%
% The following configuration options affect the behavior of
% FT_STATISTICS_MONTECARLO which is called to do the first level
% permutations for each subject. Individual subject results are returned in
% stat.ustat indexed according to stat.unqUvar.
%   cfg.precondition     = {'no', 'before', 'after'} 
%   cfg.correctm         = string, apply multiple-comparison correction, 'no', 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
%   cfg.tail             = number, -1, 1 or 0 (default = 0)
%   cfg.correcttail      = string, correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
%   cfg.wvar             = number or list with indices, within-cell variable(s)
%   cfg.cvar             = number or list with indices, control variable(s)
%   cfg.feedback         = string, 'gui', 'text', 'textbar' or 'no' (default = 'text')
%   cfg.randomseed       = string, 'yes', 'no' or a number (default = 'yes')
%
% The statistic that is computed for each sample in each random reshuffling
% of the data is specified as
%   cfg.statistic       = 'gcmi'  Gaussian-Copula Mutual Information
%  NOTE: in fact, any statistic which returns a positive-number would be fine, right?
% with extra options in 
%   cfg.gcmi (see ft_statfun_gcmi)
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS, FT_SOURCESTATISTICS,
%          FT_STATISTICS_MONTECARLO
%
% Copyright (C) 2017, Robin Ince
% Copyright (C) 2005-2015, Robert Oostenveld
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

ft_hastoolbox('prevalence-permutation',1);


ft_preamble randomseed; % deal with the user specified random seed

cfg = ft_checkconfig(cfg, 'required', {'statistic', 'ivar', 'uvar', 'numrandomization1', 'numrandomization2'});

cfg.precondition = ft_getopt(cfg, 'precondition', []);
cfg.alpha        = ft_getopt(cfg, 'alpha',      0.05);

% for backward compatibility
if size(design,2)~=size(dat,2)
  design = transpose(design);
end

% call ft_statistics_montecarlo to do permutations for each subject (uvar)
unqU = unique(cfg.uvar);
NunqU = numel(unqU);


tmpcfg           = removefields(cfg, {'uvar' 'numrandomization1' 'numrandomization2'});

% ensure the cfg to ft_statistics_montecarlo to have the correct options
tmpcfg.keeprandomizations = 'yes';
tmpcfg.resampling = 'permutation';
tmpcfg.numrandomization = cfg.numrandomization1;

% save single subject statistics structures to avoid needing to repeat it
stat = [];
stat.unqUvar = unqU;
stat.ustat = cell(1,NunqU);
% build permutations matrix

Nrand = cfg.numrandomization1;
perms = zeros(size(dat,1),NunqU,Nrand+1);
for ui=1:NunqU
  idx   = find(design(cfg.uvar,:)==unqU(ui));
  ustat = ft_statistics_montecarlo(tmpcfg, dat(:,idx), design(:,idx));
  perms(:,ui,1)     = ustat.stat;
  perms(:,ui,2:end) = ustat.statrand;
  % remove perms from structure returned for each subject
  rmfield(ustat,'statrand');
  stat.ustat{ui} = ustat;
end

% pass to prevalenceCore
[results, params] = prevalenceCore(perms, cfg.numrandomization2, cfg.alpha);


% results:      per-voxel analysis results
%   .puGN         uncorrected p-values for global null hypothesis         (Eq. 24)
%   .pcGN         corrected p-values for global null hypothesis           (Eq. 26)
%   .puMN         uncorrected p-values for majority null hypothesis       (Eq. 19)
%   .pcMN         corrected p-values for majority null hypothesis         (Eq. 21)
%   .gamma0u      uncorrected prevalence lower bounds                     (Eq. 20)
%   .gamma0c      corrected prevalence lower bounds                       (Eq. 23)
%   .aTypical     median values of test statistic where pcMN <= alpha     (Fig. 4b)

% The 'majority null hypothesis' referenced here is a special case of the
% prevalence null hypothesis (Eq. 17), where the critical value is gamma0 =
% 0.5. It describes the situation where there is no effect in the majority
% of subjects in the population. Rejecting it allows to infer that there is
% an effect in the majority of subjects in the population. aTypical is only
% defined where the (spatially extended) majority null hypothesis can be
% rejected. Compare Fig. 4b and 'What does it mean for an effect to be
% "typical" in the population?' in the Discussion of Allefeld, Goergen and
% Haynes (2016).

% median value of effect size where majority show an effect
% (prevalance > 0.5)
stat.stat = results.aTypical;
stat.mask = isfinite(stat.stat); % already nan masked
% corrected p-values for majority null hypothesis

stat.prob       = results.pcMN;
stat.probglobal = results.pcGN;
stat.prevbound  = results.gamma0c;

ft_postamble randomseed; % deal with the potential user specified randomseed



