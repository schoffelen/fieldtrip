function test_curry

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_filetype ft_read_header ft_read_data ft_read_event
% DATA private

%%
% this is where the data is installed on the DCCN compute cluster
% and where the nightly regression testing is done

cd(dccnpath('/project/3031000.02/test/original/curry'));

%%

filename = {
  './Epoched/Viscpt Epochs.cdt.dpa'
  './Epoched/Viscpt Epochs.cdt'
  './MEGandEEG/MEG Data.cdt'
  './MEGandEEG/MEG Data.cdt.dpa'
  './MEG/MEG Oddball Data.cdt.cef'
  './MEG/MEG Oddball Data.cdt.dpa'
  './MEG/MEG Oddball Data.cdt'
  './EEG/Curry 8 Data.cdt.dpa'
  './EEG/Curry 8 Data.cdt.cef'
  './EEG/Curry 7 Data.cef'
  './EEG/Curry 8 Data.cdt'
  './EEG/Curry 7 Data.dat'
  %  './EEG/Curry 7 Data.rs3' % this one is not detected/specified in ft_filetype
  './EEG/Curry 7 Data.dap'
  };


for i=1:numel(filename)
  assert(~isequal(ft_filetype(filename{i}), 'unknown'), 'file format not detected properly');
end

%%

filename = {
  './Epoched/Viscpt Epochs.cdt'
  './MEGandEEG/MEG Data.cdt'
  './MEG/MEG Oddball Data.cdt'
  './EEG/Curry 8 Data.cdt'
  './EEG/Curry 7 Data.dat'
  };

for i=1:numel(filename)
  disp(['==========' filename{i} '==========']);
  hdr = ft_read_header(filename{i});
  disp('%% header');
  disp(hdr);
  dat = ft_read_data(filename{i});
  disp('%% data');
  disp(size(dat));
  evt = ft_read_event(filename{i});
  disp('%% events');
  disp(evt);
end


%%

% these are files that were donated with https://github.com/fieldtrip/fieldtrip/issues/2446

filename = {
    './issue2446/EEGOnly.cdt'
  % './issue2446/EEGOnly.cdt.dpa'
    './issue2446/MEG + EEG + Oth.cdt'
  % './issue2446/MEG + EEG + Oth.cdt.cef'
  % './issue2446/MEG + EEG + Oth.cdt.dpa'
    './issue2446/MEG + Oth.cdt'
  % './issue2446/MEG + Oth.cdt.cef'
  % './issue2446/MEG + Oth.cdt.dpa'
    './issue2446/MEG_A + MEG_B + Oth.cdt'
  % './issue2446/MEG_A + MEG_B + Oth.cdt.ceo'
  % './issue2446/MEG_A + MEG_B + Oth.cdt.dpo'
  };

for i=1:numel(filename)
  disp(['==========' filename{i} '==========']);
  hdr = ft_read_header(filename{i});
  disp('%% header');
  disp(hdr);
  dat = ft_read_data(filename{i});
  disp('%% data');
  disp(size(dat));
  evt = ft_read_event(filename{i});
  disp('%% events');
  disp(evt);
end
