function [actpas_cds,OOR_cds] = loadsave_oor_cds(monkey,datecode,do_save)

%% Take either 3 cds files or datecode
if(~ischar(datecode))
    error('Input type error: expected datecode string')
end

%% Load CDS files
% set up input structs
root_folder = ['C:\Users\rhc307\Projects\limblab\data-raeed\ForceKin\' monkey filesep datecode filesep];
fname_prefix = [monkey '_' datecode];

lab=6;
ranBy='ranByRaeed';
monkey=['monkey' monkey];
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
actpas_cds.loadOpenSimData(folder);

% OOR_cds = commonDataStructure();
% OOR_cds.file2cds([folder filesep fname_prefix OOR_fname],      ranBy,array,monkey,lab,'ignoreJumps','taskOOR',mapfile);
% OOR_cds.file2cds([folder filesep fname_prefix extra_OOR_fname],ranBy,array,monkey,lab,'ignoreJumps','taskOOR',mapfile);

%% Save files
if(do_save)
    save_folder = [root_folder 'CDS\'];

    if(~isdir(save_folder))
        mkdir(save_folder);
    end

    save([save_folder fname_prefix '_CDS.mat'],'actpas_cds','OOR_cds','-v7.3');
end

