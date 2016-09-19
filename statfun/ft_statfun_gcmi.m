function [stat, cfg, dat] = ft_statfun_gcmi(cfg, dat, design)

% FT_STATFUN_GCMI computes mutual information between the dependent variable
% and a discrete-valued design vector.
%
% configuration-options
%  cfg.preconditionflag = 0 (default), or 1, performs Gaussian copula transform
%    Preconditioning is computationally efficient, because for given data
%    it needs to be done only once.

cfg.mi               = ft_getopt(cfg, 'mi', []);
cfg.preconditionflag = ft_getopt(cfg, 'preconditionflag', false);

% get this one as a local variable only
tra                  = ft_getopt(cfg, 'tra', []);

% check the validity of the design
if ~all(ismember(design(cfg.ivar,:), 1:max(design(cfg.ivar,:)))),
  error('the design vector is ill-specified');
end
N      = max(design(cfg.ivar,:));
Y      = design(cfg.ivar, :)' - 1; % convention of the lower level code

if cfg.preconditionflag,
  % deal with planar gradient MEG data
  if isfield(cfg, 'channel') && ft_senstype(cfg.channel, 'meg_planar')
    cfg.mi.combineplanar = ft_getopt(cfg.mi, 'combineplanar', 'no');
    
    if istrue(cfg.mi.combineplanar)
      % find the combination of horizontal and vertical channels that should be combined
      planar      = ft_senslabel(ft_senstype(cfg.channel), 'output', 'planarcombined');
      
      [~, sel_dH] = match_str(planar(:,1), cfg.channel);  % indices of the horizontal channels
      [~, sel_dV] = match_str(planar(:,2), cfg.channel);  % indices of the vertical   channels
      
      % find the other channels that are present in the data
      sel_other = setdiff(1:length(cfg.channel), [sel_dH(:)' sel_dV(:)']);
      lab_other = cfg.channel(sel_other);
      
      
      if length(sel_dH)~=length(sel_dV)
        error('not all planar channel combinations are complete')
      end
      
      % FIXME check for completeness, and that only planar channels are
      % in input
      
      % define the channel names after combining the planar combinations
      % they should be sorted according to the order of the planar channels in the data
      [~, sel_planar] = match_str(cfg.channel(sel_dH),planar(:,1));
      lab_comb        = planar(sel_planar,end);
      
      tok     = tokenize(cfg.dimord, '_');
      chandim = find(strcmp(tok, 'chan')); % FIXME works only with chandim==1
      
      indx_in  = reshape((1:size(dat,1)), cfg.dim);
      newdim   = cfg.dim;
      newdim(chandim) = numel(lab_comb)+numel(lab_other);
      indx_out = reshape((1:prod(newdim)), newdim);
      
      % create a prod(cfg.dim) x prod(newdim) sparse matrix that instructs how to combine rows of dat
      x = zeros(0,1);
      y = zeros(0,1);
      for k = 1:numel(lab_comb)
        x = cat(1,x,reshape(cat(1,indx_in(sel_dH(k),:),indx_in(sel_dV(k),:)),[],1));
        y = cat(1,y,reshape(cat(1,indx_out(k,:),indx_out(k,:)),[],1));
      end
      for k = numel(lab_comb)+(1:numel(lab_other))
        x = cat(1,x,reshape(indx_in(sel_other(k-numel(lab_comb)),:),[],1));
        y = cat(1,y,reshape(indx_out(k,:),[],1));
      end
      tra = sparse(x,y,ones(numel(x,1)));
      
    end
  end
  
  if istrue(cfg.mi.addgradient)
    if isempty(tra)
      tra = speye(size(dat,1));
    end
    
    tok     = tokenize(cfg.dimord, '_');
    timedim = find(strcmp(tok, 'time')); % FIXME works only with chandim==1
    if timedim~=2,
      error('cfg.mi.addgradient is not supported for the given data');
    end
    
    fprintf('computing the temporal derivative\n');
    dat_gradient = zeros(size(dat));
    for k = 1:size(dat,2)
      dat_gradient(:,k) = reshape(ft_preproc_derivative(reshape(dat(:,k),cfg.dim)),[],1);
    end
    dat = cat(1,dat,dat_gradient);
    clear dat_gradient;
    tra = cat(1,tra,tra);
  end
    
  if istrue(cfg.mi.combineplanar)
    % these can now be updated (i.e. only after the addgradient part)
    cfg.dim     = newdim;
    cfg.channel = [lab_comb(:);lab_other(:)];
  end
  cfg.tra     = tra; % give it back to the cfg, for re-use

  
  fprintf('performing the copula-transform\n');
  % FIXME here we should deal with NaNs in the data, note that these can be different across rows (althought that is rare), which prevents the possibility of doing the transform in a single call.
  dat  = copnorm(dat')';
  stat = [];
end

if isempty(tra)
  tra = speye(size(dat,1));
end

% convert the design to a boolean matrix
Y = indexed2boolean(Y);

% remove the class-specific mean only once (is not done anymore in mi_gd2
[dat, class_means] = centercolumns(dat,Y);


Nvar     = full(sum(tra,1));
uNvar    = unique(Nvar);
for k = 1:numel(uNvar)
  sel_col  = Nvar==uNvar(k);
  sel_row  = full(sum(tra(:,sel_col),2)>0);

  mi = mi_gd_vectorized(dat(sel_row,:)',Y,tra(sel_row,sel_col), true, true, class_means(sel_row,:)');
  stat.stat(sel_col,1) = mi;
end

% mi = zeros(size(tra,2),1);
% for k = 1:size(tra,2)
%   mi(k) = mi_gd2(dat(tra(:,k)>0,:)',Y, N, true, true, class_means(tra(:,k)>0,:)'); %JM's version that works with persistent variables, to speed up some computations
% end
% stat.stat = mi;

if nargout>2,
  % add the class-specific means back to the data
  dat = addclassmeans(dat,Y,class_means);
end

function [dat, class_means] = centercolumns(dat, design)

%udesign     = unique(design);
%class_means = zeros(size(dat,1),numel(udesign));
class_means = zeros(size(dat,1),size(design,2));
%for k = 1:numel(udesign)
for k = 1:size(design,2)
  %sel = design==udesign(k);
  sel = design(:,k);
  class_means(:,k) = mean(dat(:,sel),2);
  dat(:,sel) = bsxfun(@minus,dat(:,sel),class_means(:,k));
end

function dat = addclassmeans(dat,design,class_means)

%udesign = unique(design);
%for k = 1:numel(udesign)
for k = 1:size(design,2)
  %sel = design==udesign(k);
  sel = design(:,k);
  dat(:,sel) = bsxfun(@plus,dat(:,sel),class_means(:,k));
end

function Y = indexed2boolean(X)

uX = unique(X);
Y  = false(numel(X),numel(uX));
for k = 1:size(Y,2)
  Y(X==uX(k),k) = true;
end
