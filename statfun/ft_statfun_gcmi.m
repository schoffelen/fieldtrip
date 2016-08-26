function [stat, cfg, dat] = statfun_anova2x2rm(cfg, dat, design)

% FT_STATFUN_GCMI computes mutual information between the dependent variable
% and a discrete-valued design vector.
%
% configuration-options
%  cfg.preconditionflag = 0 (default), or 1, performs Gaussian copula transform
%    Preconditioning is computationally efficient, because for given data
%    it needs to be done only once.

cfg.preconditionflag = ft_getopt(cfg, 'preconditionflag', false);

% check the validity of the design
if ~all(ismember(design(cfg.ivar,:), 1:max(design(cfg.ivar,:)))),
  error('the design vector is ill-specified');
end
N      = max(design(cfg.ivar,:));
y      = design(cfg.ivar, :) - 1; % convention of the lower level code

if cfg.preconditionflag,
  fprintf('performing the copula-transform\n');
  % FIXME here we should deal with NaNs in the data, note that these can be different across rows (althought that is rare), which prevents the possibility of doing the transform in a single call.
  dat  = copnorm(dat')';
  stat = [];
end

mi = zeros(size(dat,1),1);
for k = 1:size(dat,1)
  mi(k) = mi_gd(dat(k,:)',y, N, true, true);
end
stat.stat = mi;

