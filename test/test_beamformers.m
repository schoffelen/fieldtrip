function test_beamformers

% WALLTIME 00:45:00
% MEM 6gb
% DEPENDENCY ft_sourceanalysis ft_inverse_lcmv

%% do some virtual channel stuff
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/source_coh_lft.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/source_diff.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/data_cmb.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/sourcemodel.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/hdm.mat'));

[maxval, maxcohindx] = max(source_coh_lft.avg.coh);
[maxval, maxpowindx] = max(source_diff.avg.pow);

cfg                   = [];
cfg.covariance        = 'yes';
cfg.channel           = 'MEG';
cfg.vartrllength      = 2;
cfg.covariancewindow  = 'all';
tlock                 = ft_timelockanalysis(cfg, data_cmb);

% this is old-style stuff. as of end 2020 there's a ft_virtualchannel
% function that does the virtualchannel creation
cfg             = [];
cfg.headmodel   = hdm;
cfg.sourcemodel = sourcemodel;
cfg.channel     = 'MEG';
cfg.singleshell.batchsize = 2000;
leadfield       = ft_prepare_leadfield(cfg, tlock);

cfg              = [];
cfg.method       = 'lcmv';
cfg.sourcemodel  = leadfield;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'no';
source1           = ft_sourceanalysis(cfg, tlock);

cfg.lcmv.fixedori = 'yes';
source2           = ft_sourceanalysis(cfg, tlock);

cfg.lcmv.weightnorm = 'arraygain';
cfg.lcmv.fixedori = 'no';
source3           = ft_sourceanalysis(cfg, tlock);

cfg.lcmv.weightnorm = 'arraygain';
cfg.lcmv.fixedori = 'yes';
source4           = ft_sourceanalysis(cfg, tlock);

cfg.lcmv.weightnorm = 'unitnoisegain';
cfg.lcmv.fixedori = 'no';
source5           = ft_sourceanalysis(cfg, tlock);

cfg.lcmv.weightnorm = 'unitnoisegain';
cfg.lcmv.fixedori = 'yes';
source6           = ft_sourceanalysis(cfg, tlock);

inside = find(source1.inside);
wl1 = zeros(source1.dim);
wl2 = zeros(source2.dim);
wl3 = wl1;
wl4 = wl2;
wl5 = wl1;
wl6 = wl2;
ww5 = wl5; % filter sum-of-squares (for unitnoisegain)
ww6 = wl6;
ll3 = wl3; % leadfield sum-of-squares (for arraygain)
ll4 = wl4;
for k = 1:numel(inside)
  [ix1,ix2,ix3]=ind2sub(source1.dim,inside(k));
  
  % vector versions
  wl1(ix1,ix2,ix3) = trace(source1.avg.filter{inside(k)}*leadfield.leadfield{inside(k)}); %rotated wtl=eye(2); Stimmt
  wl3(ix1,ix2,ix3) = trace(source3.avg.filter{inside(k)}*leadfield.leadfield{inside(k)}); %?
  wl5(ix1,ix2,ix3) = trace(source5.avg.filter{inside(k)}*leadfield.leadfield{inside(k)}); %?
  
  % scalar versions
  wl2(ix1,ix2,ix3) = source2.avg.filter{inside(k)}*leadfield.leadfield{inside(k)}*source2.avg.ori{inside(k)}; %wtl=1; Stimmt
  wl4(ix1,ix2,ix3) = source4.avg.filter{inside(k)}*leadfield.leadfield{inside(k)}*source4.avg.ori{inside(k)}; %wtl=sqrt(ltl); Stimmt
  wl6(ix1,ix2,ix3) = source6.avg.filter{inside(k)}*leadfield.leadfield{inside(k)}*source6.avg.ori{inside(k)}; %Stimmt in line with wtw=1 
  
  % filter norm, as per checking the unitnoisegain constraint
  ww5(ix1,ix2,ix3) = trace(source5.avg.filter{inside(k)}*source5.avg.filter{inside(k)}'); %rotated wtw=eye(3); Stimmt NICHT
  ww6(ix1,ix2,ix3) = source6.avg.filter{inside(k)}*source6.avg.filter{inside(k)}'; %wtw=1; Stimmt
  
  % leadfield norm, as per checking the arraygain constraint
  ll3(ix1,ix2,ix3) = sum(sum(leadfield.leadfield{inside(k)}.^2));
  ll4(ix1,ix2,ix3) = sum((leadfield.leadfield{inside(k)}*source4.avg.ori{inside(k)}).^2);
end
