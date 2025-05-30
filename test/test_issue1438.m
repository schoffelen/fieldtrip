function test_issue1438

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY read_edf
% DATA private

filename = dccnpath('/project/3031000.02/test/issue1438/test-saline-run3.edf');
hdr = ft_read_header(filename);
event = ft_read_event(filename, 'header', hdr, 'detectflank', 'up', 'chanindx', [6:13], 'threshold', '5*nanmedian');

assert(numel(event)==587);