%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pdTable = get_targ_PDs(trial_data,params)
%
%   Gets PD table for given out_signal. You need to define the out_signal
% and move_corr parameters at input.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .out_signals  : which signals to calculate PDs for
%       .out_signal_names : cell array of names of signals to be used as signalID pdTable
%                           default - empty
%       .trial_idx    : trials to use.
%                         DEFAULT: [1,length(trial_data]
%       .conditions   : which conditions to calculate PDs on
%                           note: must be directional conditions
%                           default - {'target_direction';'force_direction'}
%       .num_boots    : # bootstrap iterations to use
%       .distribution : distribution to use. See fitglm for options
%
% OUTPUTS:
%   pdTable : calculated PD table with CIs
%
% Written by Raeed Chowdhury. Updated Nov 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pdTable = get_targ_PDs(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
out_signals      =  [];
out_signal_names = {};
trial_idx        =  1:length(trial_data);
conditions      = {'target_direction';'force_direction'};
num_boots        =  1000;
distribution = 'Normal';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented parameters
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
trial_data = trial_data(trial_idx);

if ~strcmpi(trial_data(1).task,'OOR')
    warning('Only tested for OOR task, behavior may be different for other tasks')
end

if isempty(out_signals), error('Need to provide output signal'); end

out_signals = check_signals(trial_data(1),out_signals);
response_var = get_vars(trial_data,out_signals);

input_var = [];
for i = 1:size(conditions,1)
    input_cond = cat(1,trial_data.(conditions{i}));
    input_var = [input_var cos(input_cond) sin(input_cond)];
end

if numel(unique(cat(1,{trial_data.monkey}))) > 1
    error('More than one monkey in trial data')
end
monkey = repmat({trial_data(1).monkey},size(response_var,2),1);
if numel(unique(cat(1,{trial_data.date}))) > 1
    date = cell(size(response_var,2),1);
    warning('More than one date in trial data')
else
    date = repmat({trial_data(1).date},size(response_var,2),1);
end
if numel(unique(cat(1,{trial_data.task}))) > 1
    task = cell(size(response_var,2),1);
    warning('More than one task in trial data')
else
    task = repmat({trial_data(1).task},size(response_var,2),1);
end

out_signal_names = reshape(out_signal_names,size(response_var,2),[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preallocate final table
dirArr = zeros(size(response_var,2),1);
dirCIArr = zeros(size(response_var,2),2);
moddepthArr = zeros(size(response_var,2),1);
moddepthCIArr = zeros(size(response_var,2),2);
pdTable = table(monkey,date,task,out_signal_names,'VariableNames',{'monkey','date','task','signalID'});
for in_signal_idx = 1:size(conditions,1)
    tab_append = table(dirArr,dirCIArr,moddepthArr,moddepthCIArr,...
                        'VariableNames',{[conditions{in_signal_idx,1} 'PD'],[conditions{in_signal_idx,1} 'PDCI'],[conditions{in_signal_idx,1} 'Moddepth'],[conditions{in_signal_idx,1} 'ModdepthCI']});
    pdTable = [pdTable tab_append];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get mean firing rates for each trial
fr = zeros(length(trial_data),size(trial_data(1).S1_spikes,2));
for trial_idx = 1:length(trial_data)
    fr(trial_idx,:) = mean(get_vars(trial_data(trial_idx),out_signals));
end

% fit tuning curves
if ~strcmpi(distribution,'normal')
    warning('GLM not implemented, using LM instead')
end
bootfunc = @(data) fitlm(data(:,2:end),data(:,1));
tic;
for uid = 1:size(fr,2)
    disp(['  Bootstrapping GLM PD computation(ET=',num2str(toc),'s).'])
    data_arr = [fr(:,uid) input_var];
    boot_tuning = bootstrp(num_boots,@(data) {bootfunc(data)}, data_arr);
    boot_coef = cell2mat(cellfun(@(x) x.Coefficients.Estimate',boot_tuning,'uniformoutput',false));

    if size(boot_coef,2) ~= 1+size(conditions,1)*2
        error('get_targ_PDs:moveCorrProblem','GLM doesn''t have correct number of inputs')
    end

    for in_signal_idx = 1:size(conditions,1)
        move_corr = conditions{in_signal_idx,1};

        dirs = atan2(boot_coef(:,1+in_signal_idx*2),boot_coef(:,in_signal_idx*2));
        %handle wrap around problems:
        centeredDirs=minusPi2Pi(dirs-circ_mean(dirs));

        pdTable.([move_corr 'PD'])(uid,:)=circ_mean(dirs);
        pdTable.([move_corr 'PDCI'])(uid,:)=prctile(centeredDirs,[2.5 97.5])+circ_mean(dirs);

        if strcmpi(distribution,'normal')
            % get moddepth
            moddepths = sqrt(sum(boot_coef(:,(2*in_signal_idx):(2*in_signal_idx+1)).^2,2));
            pdTable.([move_corr 'Moddepth'])(uid,:)= mean(moddepths);
            pdTable.([move_corr 'ModdepthCI'])(uid,:)= prctile(moddepths,[2.5 97.5]);
        else
            pdTable.([move_corr 'Moddepth'])(uid,:)= -1;
            pdTable.([move_corr 'ModdepthCI'])(uid,:)= [-1 -1];
        end
    end
end
