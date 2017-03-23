function [trial_data_actpas,trial_data_OOR] = loadsave_oor_td(datecode,do_save)

%% Take either 3 cds files or datecode
if(~ischar(datecode))
    error('Input type error: expected datecode string')
end

%% Load CDS files
% set up input structs
root_folder = ['C:\Users\rhc307\Projects\limblab\data-raeed\ForceKin\Han\' datecode filesep];
fname_prefix = ['Han_' datecode];

lab=6;
ranBy='ranByRaeed';
monkey='monkeyHan';
task='taskOOR';
array='arrayLeftS1Area2';
folder=[root_folder 'preCDS'];
mapfile='mapFileC:\Users\rhc307\Projects\limblab\data-raeed\ForceKin\OutOutReach\Han\mapfile\left_S1\SN 6251-001459.cmp';

actpas_fname =  '_COactpas_area2EMG_001';
OOR_fname =     '_OOR_25N_area2EMG_002';

extra_actpas_fname =    '_COactpas_EMGextra_001';
extra_OOR_fname =       '_OOR_25N_EMGextra_002';

actpas_cds = commonDataStructure();
actpas_cds.file2cds([folder filesep fname_prefix actpas_fname],      ranBy,array,monkey,lab,'ignoreJumps','taskCObump',mapfile);
actpas_cds.file2cds([folder filesep fname_prefix extra_actpas_fname],ranBy,array,monkey,lab,'ignoreJumps','taskCObump',mapfile);

OOR_cds = commonDataStructure();
OOR_cds.file2cds([folder filesep fname_prefix OOR_fname],      ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);
OOR_cds.file2cds([folder filesep fname_prefix extra_OOR_fname],ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);

%% Extract trial data from cds files
params.array_alias = {'LeftS1Area2','S1'};
% params.exclude_units = [255];
params.event_list = {'ctrHoldBump';'bumpTime';'bumpDir'};
params.trial_results = {'R','A','F','I'};
meta = struct('task','COactpas');
params.meta = meta;
trial_data_actpas = parseFileByTrial(actpas_cds,params);

params.event_list = {'forceDir'};
params.meta.task = 'OOR';
trial_data_OOR = parseFileByTrial(OOR_cds,params);

% trial_data = cat(2,trial_data_actpas,trial_data_OOR);


%% Save files
if(do_save)
    save_folder = [root_folder 'CDS\'];

    if(~isdir(save_folder))
        mkdir(save_folder);
    end

    save([save_folder fname_prefix '_CDS_noLFP.mat'],'actpas_cds','OOR_cds','-v7.3');

    save_folder = [root_folder 'TD\'];

    if(~isdir(save_folder))
        mkdir(save_folder);
    end

    save([save_folder fname_prefix '_TD.mat'],'trial_data_actpas','trial_data_OOR','-v7.3');
end

