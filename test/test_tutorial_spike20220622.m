function test_tutorial_spike20220622

% WALLTIME 00:30:00
% MEM 6gb
% DEPENDENCY ft_read_spike ft_spike_select ft_spike_waveform ft_spike_maketrials ft_spike_plot_jpsth
% DATA public

cd(dccnpath('/project/3031000.02/external/download/tutorial/spike/'));

%%

spike = ft_read_spike('p029_sort_final_01.nex');

cfg = [];
cfg.spikechannel = {'sig002a_wf', 'sig003a_wf'}; % select only the two single-units
spike = ft_spike_select(cfg, spike);

%% ## Computing average waveforms

cfg             = [];
cfg.fsample     = 40000;
cfg.interpolate = 1; % keep the density of samples as is
[wave, spikeCleaned] = ft_spike_waveform(cfg,spike);

% The resulting wave structure has the following content
for k = [1 2]
  figure,
  sl = squeeze(wave.dof(k,:,:))>1000; % only keep samples with enough spikes
  plot(wave.time(sl), squeeze(wave.avg(k,:,sl)),'k') % factor 10^6 to get microseconds
  hold on

  % plot the standard deviation
  plot(wave.time(sl), squeeze(wave.avg(k,:,sl))+sqrt(squeeze(wave.var(k,:,sl))),'k--')
  plot(wave.time(sl), squeeze(wave.avg(k,:,sl))-sqrt(squeeze(wave.var(k,:,sl))),'k--')

  axis tight
  set(gca,'Box', 'off')
  xlabel('time')
  ylabel('normalized voltage')
end

%% ## Adding trigger event information to spike structure

event = ft_read_event('p029_sort_final_01.nex');

cfg          = [];
cfg.dataset  = 'p029_sort_final_01.nex';
cfg.trialfun = 'trialfun_stimon';
cfg = ft_definetrial(cfg);

cfg.timestampspersecond =  spike.hdr.FileHeader.Frequency; % 40000
spikeTrials = ft_spike_maketrials(cfg,spike);

plot(spikeTrials.timestamp{1}, spikeTrials.time{1}, '.')
axis([1e6 5.5e6 -3 5])

spikeTrials.trialtime(1:5,:)

spikeTrials.cfg.trl(1:5,:)

cfg                     = [];
hdr                     = ft_read_header('p029_sort_final_01.nex');
cfg.trl                 = [0 hdr.nSamples*hdr.TimeStampPerSample 0];
cfg.timestampspersecond =  spike.hdr.FileHeader.Frequency; % 40000
spike_notrials   = ft_spike_maketrials(cfg,spike);

dat = ft_checkdata(spikeTrials,'datatype', 'raw', 'fsample', 1000)

dat.trial{1}(:,4000:4004)
spike_converted = ft_checkdata(dat,'datatype', 'spike')


%% ## Characterizing inter-spike-interval (ISI) distributions

cfg       = [];
cfg.bins  = [0:0.0005:0.1]; % use bins of 0.5 milliseconds
cfg.param = 'coeffvar'; % compute the coefficient of variation (sd/mn of isis)
isih = ft_spike_isi(cfg,spikeTrials);

for k = [1 2] % only do for the single-units
  cfg              = [];
  cfg.spikechannel = isih.label{k};
  cfg.interpolate  = 5; % interpolate at 5 times the original density
  cfg.window       = 'gausswin'; % use a gaussian window to smooth
  cfg.winlen       = 0.004; % the window by which we smooth has size 4 by 4 ms
  cfg.colormap     = jet(300); % colormap
  cfg.scatter      = 'no'; % do not plot the individual isis per spike as scatters
  figure, ft_spike_plot_isireturn(cfg,isih)
end

filename    = 'tt6_7.t'
cfg         = [];
cfg.dataset = filename;
spike2 = ft_read_spike(cfg.dataset)

% convert timestamps to seconds
cfg                     = [];
cfg.trl                 = [0 max(spike2.timestamp{1})+1 0];
cfg.timestampspersecond = 10^6;
spike2Trial = ft_spike_maketrials(cfg,spike2);

% run the isi histogram
cfg      = [];
cfg.bins = [0:0.001:0.2];
isih = ft_spike_isi(cfg,spike2Trial);

% plot the isi histogram with the Poincare return map
cfg             = [];
cfg.interpolate = 5;
cfg.window      = 'gausswin';
cfg.winlen      = 0.005;
cfg.scatter     = 'no';
cfg.colormap    = jet(300);
figure, ft_spike_plot_isireturn(cfg,isih)

%% ## Computing spike densities and peri-stimulus time histograms (PSTHs)

cfg             = [];
cfg.binsize     =  0.1; % if cfgPsth.binsize = 'scott' or 'sqrt', we estimate the optimal bin size from the data itself
cfg.outputunit  = 'rate'; % give as an output the firing rate
cfg.latency     = [-1 3]; % between -1 and 3 sec.
cfg.vartriallen = 'yes'; % variable trial lengths are accepted
cfg.keeptrials  = 'yes'; % keep the psth per trial in the output
psth = ft_spike_psth(cfg,spikeTrials);

cfg         = [];
cfg.binsize =  [0.05];
cfg.latency = [-1 3];
psth = ft_spike_psth(cfg,spikeTrials);

cfg              = [];
cfg.topplotfunc  = 'line'; % plot as a line
cfg.spikechannel = spikeTrials.label([1 2]);
cfg.latency      = [-1 3];
cfg.errorbars    = 'std'; % plot with the standard deviation
cfg.interactive  = 'no'; % toggle off interactive mode
figure, ft_spike_plot_raster(cfg,spikeTrials, psth)


cfg         = [];
cfg.latency = [-1 3];
cfg.timwin  = [-0.025 0.025];
cfg.fsample = 1000; % sample at 1000 hz
sdf = ft_spikedensity(cfg,spikeTrials);

cfg              = [];
cfg.topplotfunc  = 'line'; % plot as a line plot
cfg.spikechannel = spikeTrials.label([1 2]); % can also select one unit here
cfg.latency      = [-1 3];
cfg.errorbars    = 'std'; % plot with standard deviation
cfg.interactive  = 'no'; % toggle off interactive mode
figure, ft_spike_plot_raster(cfg,spikeTrials, sdf)

cfg         = [];
cfg.latency = [-1 3];
cfg.timwin  = [-0.1 0.1];
cfg.fsample = 1000; % sample at 1000 hz
[sdf, sdfdata] = ft_spikedensity(cfg,spikeTrials);

sdfdata.trial{1}(:,3000:3005)

%% ## Computing average firing rates and noise correlations

cfg            = [];
cfg.latency    = [0.3 max(spikeTrials.trialtime(:))]; % sustained response period
cfg.keeptrials = 'yes';
rate = ft_spike_rate(cfg,spikeTrials);

[R,P] = corrcoef(rate.trial)

R % Pearson's R

P % probability

%% ## Computing cross-correlations between spike trains

spike = ft_read_spike('p029_sort_final_01.nex');

cfg              = [];
cfg.spikechannel = {'sig001U_wf', 'sig002U_wf', 'sig003U_wf', 'sig004U_wf'}; % select only the MUA
spike = ft_spike_select(cfg, spike);

cfg          = [];
cfg.dataset  = 'p029_sort_final_01.nex';
cfg.trialfun = 'trialfun_stimon';
cfg = ft_definetrial(cfg);

cfg.timestampspersecond =  spike.hdr.FileHeader.Frequency; % 40000
spikeTrials = ft_spike_maketrials(cfg,spike);

% and then compute the cross-correlogram (and the shift-predictor cross-correlogram) by
%
cfg             = [];
cfg.maxlag      = 0.2; % maximum 200 ms
cfg.binsize     = 0.001; % bins of 1 ms
cfg.outputunit  = 'proportion'; % make unit area
cfg.latency     = [-2.5 0];
cfg.vartriallen = 'no'; % do not allow variable trial lengths
cfg.method      = 'xcorr'; % compute the normal cross-correlogram
Xc = ft_spike_xcorr(cfg,spikeTrials);

% compute the shuffled correlogram
cfg.method      = 'shiftpredictor'; % compute the shift predictor
Xshuff = ft_spike_xcorr(cfg,spikeTrials);

iCmb = 3;
jCmb = 4;
figure
xcSmoothed = conv(squeeze(Xc.xcorr(iCmb,jCmb,:)),ones(1,5)./5,'same'); % do some smoothing
hd = plot(Xc.time(3:end-2),xcSmoothed(3:end-2),'k'); % leave out borders (because of smoothing)
hold on
xcSmoothed = conv(squeeze(Xshuff.shiftpredictor(iCmb,jCmb,:)),ones(1,5)./5,'same');
plot(Xc.time(3:end-2),xcSmoothed(3:end-2),'r')
hold on
xlabel('delay')
ylabel('proportion of coincidences')
title([Xc.label{iCmb} Xc.label{jCmb}])
axis tight

%% ## The joint peri stimulus time histogram

cfg         = [];
cfg.latency = [-1 3];
cfg.timwin  = [-0.025 0.025];
cfg.fsample = 200;
[sdf] = ft_spikedensity(cfg,spikeTrials);

% compute the joint psth
cfg               = [];
cfg.normalization = 'no';
cfg.channelcmb    = spikeTrials.label(3:4);
cfg.method        = 'jpsth';
jpsth = ft_spike_jpsth(cfg,sdf);

cfg.method = 'shiftpredictor';
jpsthShuff = ft_spike_jpsth(cfg,sdf);

% subtract the predictor
jpsthSubtr = jpsth;
jpsthSubtr.jpsth = jpsth.jpsth-jpsthShuff.shiftpredictor;

cfg        = [];
figure
ft_spike_plot_jpsth(cfg,jpsth)
figure, ft_spike_plot_jpsth(cfg,jpsthShuff)
figure, ft_spike_plot_jpsth(cfg,jpsthSubtr)

