function test_issue2012

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_virtualchannel
% DATA private

filename = dccnpath('/project/3031000.02/test/issue2012.mat');
load(filename);

cfg.pos = cfg.pos(1:100,:); % otherwise too memory greedy, not needed for a proof of principle
v=ft_virtualchannel(cfg,TFData,source_av);