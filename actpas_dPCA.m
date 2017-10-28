%% setup
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

% get still handle data (no word for start of center hold)
minCH = min(cat(1,td.ctrHold));
bin_size = td(1).bin_size;
still_bins = floor(minCH/bin_size);

% Get td_act and td_pas
num_bins_before = floor(still_bins/2);
num_bins_after = 30;

% prep td_act
[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_movement_on',-num_bins_before},{'idx_movement_on',num_bins_after});
% clean nans out...?
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);

% prep td_pas
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',-num_bins_before},{'idx_bumpTime',num_bins_after});
% move bumpDir into target_direction for passive td
if floor(td_pas(1).bumpDir) == td_pas(1).bumpDir
    % probably in degrees
    multiplier = pi/180;
else
    warning('bumpDir may be in radians')
    multiplier = 1;
end
for trial = 1:length(td_pas)
    td_pas(trial).target_direction = td_pas(trial).bumpDir*multiplier;
end

% concatenate tds
td = cat(2,td_act,td_pas);

%% try dPCA
td = sqrtTransform(td,'S1_spikes');
td = smoothSignals(td,struct('signals',{{'S1_spikes'}}));
[td,dPCA_info] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'S1_spikes'}},'num_dims',8,'do_plot',true,...
                                                                    'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                    'dpca_plot_fcn',@dpca_plot_actpas));

%% Try fabricating trial_data with linear models based on handle kinematics and force
% get models for force and velocity from actpas data
[td,model_info] = getModel(td,struct('model_type','linmodel',...
    'model_name','S1_handle','in_signals',{{'vel';'force'}},...
    'out_signals',{'S1_spikes'}));

% [td,dPCA_info] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'linmodel_S1_handle'}},'num_dims',4,'do_plot',true,...
%                                                                     'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
%                                                                     'dpca_plot_fcn',@dpca_plot_actpas));

%% Try fabricating trial_data with linear models based on muscles
% get models for force and velocity from actpas data
% opensim_idx = find(contains(td(1).opensim_names,'_moment'));
% opensim_idx = find(contains(td(1).opensim_names,'_vel'));
opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
% opensim_idx = find(contains(td(1).opensim_names,'_vel') & ~contains(td(1).opensim_names,'wrist') & ~contains(td(1).opensim_names,'radial'));
% opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
% opensim_idx = find(contains(td(1).opensim_names,'_vel') | contains(td(1).opensim_names,'_moment'));
% opensim_idx = find( (contains(td(1).opensim_names,'_vel') | contains(td(1).opensim_names,'_moment')) & ~contains(td(1).opensim_names,'wrist') & ~contains(td(1).opensim_names,'radial'));
[td,model_info] = getModel(td,struct('model_type','linmodel',...
    'model_name','S1_muscle','in_signals',...
    {{'opensim',opensim_idx}},...
    'out_signals',{'S1_spikes'}));

td = getPCA(td,struct('signals',{{'linmodel_S1_muscle'}}));
[td,dPCA_info] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'linmodel_S1_muscle_pca',1:length(opensim_idx)}},'num_dims',8,'do_plot',true,...
                                                                    'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                    'dpca_plot_fcn',@dpca_plot_actpas));
% [td,dPCA_info] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'linmodel_S1_muscle'}},'num_dims',8,'do_plot',true,...
%                                                                     'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
%                                                                     'dpca_plot_fcn',@dpca_plot_actpas));

muscle_noise_covar = getNoiseCovar(td,struct('actual_signals',{{'S1_spikes'}},'modeled_signals',{{'linmodel_S1_muscle'}},...
                                                'do_plot',false));
td = addCorrelatedNoise(td,struct('signals',{{'linmodel_S1_muscle'}},'noise_covar',muscle_noise_covar));
% td = smoothSignals(td,struct('signals',{{'linmodel_S1_muscle_noisy'}}));
[td,dPCA_info] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'linmodel_S1_muscle_noisy'}},'num_dims',8,'do_plot',true,...
                                                                    'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                    'dpca_plot_fcn',@dpca_plot_actpas));
