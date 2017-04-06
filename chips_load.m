%%
datecode = '20170323';
root_folder = ['C:\Users\rhc307\Projects\limblab\data-raeed\ForceKin\Chips\' datecode filesep];
fname_prefix = ['Chips_' datecode];

lab=6;
ranBy='ranByRaeed';
monkey='monkeyChips';
array='arrayLeftS1Area2';
folder=[root_folder 'preCDS'];
mapfile='mapFileC:\Users\rhc307\Projects\limblab\data-raeed\Meta\Mapfiles\Map files\Chips\left_S1\SN 6251-001455.cmp';

actpas_fname =  '_COactpas_area2_001';

actpas_cds = commonDataStructure();
actpas_cds.file2cds([folder filesep fname_prefix actpas_fname],      ranBy,array,monkey,lab,'ignoreJumps','taskCObump',mapfile);

%% Extract trial data from cds files
params.array_alias = {'LeftS1Area2','S1'};
params.exclude_units = [0,255];
params.event_list = {'ctrHoldBump';'bumpTime';'bumpDir'};
params.trial_results = {'R','A','F','I'};
meta = struct('task','COactpas');
params.meta = meta;
trial_data_actpas = parseFileByTrial(actpas_cds,params);

% trial_data = cat(2,trial_data_actpas,trial_data_OOR);