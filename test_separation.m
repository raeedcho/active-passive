%% Get PCA for act vs pas
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
td_act = binTD(td_act,15);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_pas = binTD(td_pas,15);
td = cat(2,td_act,td_pas);

%% test separability for actual spikes
td = sqrtTransform(td,'S1_spikes');
% td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
td = getPCA(td,struct('signals',{{'S1_spikes'}}));
test_sep(td,struct('signals',{{'S1_pca'}},'do_plot',true))

%% Try fabricating trial_data with linear models based on handle kinematics and force
% get models for force and velocity from actpas data
[td,model_info] = getModel(td,struct('model_type','linmodel',...
    'model_name','S1_handle','in_signals',{{'vel';'force'}},...
    'out_signals',{'S1_spikes'}));

td = getPCA(td,struct('signals',{{'linmodel_S1_handle'}}));
test_sep(td,struct('signals',{{'linmodel_S1_handle_pca',1:4}},'do_plot',true))

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
% test_sep(td,struct('signals',{{'linmodel_S1_pca',1:length(opensim_idx)}},'do_plot',true))
test_sep(td,struct('signals',{{'linmodel_S1_muscle_pca'}},'do_plot',true))

%% get boostrapped separability values
n_boot = 10;
bootsep_true = bootstrp(n_boot,@(x) test_sep(x,struct('signals',{{'S1_pca'}})),td);
bootsep_handle = bootstrp(n_boot,@(x) test_sep(x,struct('signals',{{'linmodel_S1_handle_pca'}})),td);
bootsep_muscle = bootstrp(n_boot,@(x) test_sep(x,struct('signals',{{'linmodel_S1_muscle_pca'}})),td);

bar([mean(bootsep_true) mean(bootsep_handle) mean(bootsep_muscle)])
hold on
sep_true = prctile(bootsep_true,[2.5 97.5]);
sep_handle = prctile(bootsep_handle,[2.5 97.5]);
sep_muscle = prctile(bootsep_muscle,[2.5 97.5]);
plot([1:3;1:3],[sep_true sep_handle sep_muscle],'-k','linewidth',2)

%% Try fabricating trial_data with straight up handle stuff
% [~,td] = getTDidx(trial_data,'result','R');
% 
% td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
% 
% [~,td_act] = getTDidx(td,'ctrHoldBump',false);
% td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
% td_act = binTD(td_act,15);
% [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
% td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
% td_pas = binTD(td_pas,15);
% td = cat(2,td_act,td_pas);
% 
% td = getPCA(td,struct('signals',{{'vel';'force'}}));
% test_sep(td,struct('signals',{{'velforce_pca'}},'do_plot',true))

%% Try fabricating trial_data with straight up joint stuff
% [~,td] = getTDidx(trial_data,'result','R');
% 
% td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
% 
% [~,td_act] = getTDidx(td,'ctrHoldBump',false);
% td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
% td_act = binTD(td_act,15);
% [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
% td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
% td_pas = binTD(td_pas,15);
% td = cat(2,td_act,td_pas);
% 
% % opensim_idx = find(contains(td(1).opensim_names,'_vel'));
% % opensim_idx = find(contains(td(1).opensim_names,'_vel') | contains(td(1).opensim_names,'_moment'));
% td = getPCA(td,struct('signals',{{'opensim',opensim_idx}}));
% test_sep(td,struct('signals',{{'opensim_pca'}},'do_plot',true))
