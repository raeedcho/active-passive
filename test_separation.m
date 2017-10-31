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

figure
subplot(2,3,4)
[actual_sep,actual_mdl] = test_sep(td,struct('signals',{{'S1_pca'}},'do_plot',true));
disp(['Actual separability - ' num2str(actual_sep)])

% then for directional separability/other view
signal_act = get_vars(td_act,{{'S1_pca'}});
signal_pas = get_vars(td_pas,{{'S1_pca'}});
% plot active as filled, passive as open
bump_colors = linspecer(4);
act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;
subplot(2,3,1)
hold all
scatter3(signal_act(:,1),signal_act(:,2),signal_act(:,3),50,bump_colors(act_dir_idx,:),'filled')
scatter3(signal_pas(:,1),signal_pas(:,2),signal_pas(:,3),100,bump_colors(pas_dir_idx,:),'o','linewidth',2)
ylim = get(gca,'ylim');
zlim = get(gca,'zlim');
set(gca,'box','off','tickdir','out')
axis equal
axis off

%% Try fabricating trial_data with linear models based on handle kinematics and force
% get models for force and velocity from actpas data
[td,model_info] = getModel(td,struct('model_type','linmodel',...
    'model_name','S1_handle','in_signals',{{'vel';'force'}},...
    'out_signals',{'S1_spikes'}));

td = getPCA(td,struct('signals',{{'linmodel_S1_handle'}}));

figure
subplot(2,3,5)
[velforce_sep,velforce_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_handle_pca',1:4}},'do_plot',true));
% [velforce_sep,velforce_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_handle_pca'}},'mdl',actual_mdl,'do_plot',true));
disp(['Velforce separability - ' num2str(velforce_sep)])

% then for directional separability/other view
signal_act = get_vars(td_act,{{'linmodel_S1_handle_pca'}});
signal_pas = get_vars(td_pas,{{'linmodel_S1_handle_pca'}});
% plot active as filled, passive as open
bump_colors = linspecer(4);
act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;
subplot(2,3,2)
hold all
scatter3(signal_act(:,1),signal_act(:,2),signal_act(:,3),50,bump_colors(act_dir_idx,:),'filled')
scatter3(signal_pas(:,1),signal_pas(:,2),signal_pas(:,3),100,bump_colors(pas_dir_idx,:),'o','linewidth',2)
ylim = get(gca,'ylim');
zlim = get(gca,'zlim');
set(gca,'box','off','tickdir','out')
axis equal
axis off

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

figure
subplot(2,3,6)
[muscle_sep,muscle_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_muscle_pca',1:length(opensim_idx)}},'do_plot',true));
% [muscle_sep,muscle_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_muscle_pca'}},'do_plot',true));
% [muscle_sep,muscle_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_muscle_pca'}},'mdl',actual_mdl,'do_plot',true));
disp(['Muscle separability - ' num2str(muscle_sep)])

% then for directional separability/other view
signal_act = get_vars(td_act,{{'linmodel_S1_muscle_pca'}});
signal_pas = get_vars(td_pas,{{'linmodel_S1_muscle_pca'}});
% plot active as filled, passive as open
bump_colors = linspecer(4);
act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;
subplot(2,3,3)
hold all
scatter3(signal_act(:,1),signal_act(:,2),signal_act(:,3),50,bump_colors(act_dir_idx,:),'filled')
scatter3(signal_pas(:,1),signal_pas(:,2),signal_pas(:,3),100,bump_colors(pas_dir_idx,:),'o','linewidth',2)
ylim = get(gca,'ylim');
zlim = get(gca,'zlim');
set(gca,'box','off','tickdir','out')
axis equal
axis off

%% get boostrapped separability values

% get correlated noise models
handle_noise_covar = getNoiseCovar(td,struct('actual_signals',{{'S1_spikes'}},'modeled_signals',{{'linmodel_S1_handle'}},...
                                                'do_plot',false));
td = addCorrelatedNoise(td,struct('signals',{{'linmodel_S1_handle'}},'noise_covar',handle_noise_covar));
td = getPCA(td,struct('signals',{{'linmodel_S1_handle_noisy'}}));
muscle_noise_covar = getNoiseCovar(td,struct('actual_signals',{{'S1_spikes'}},'modeled_signals',{{'linmodel_S1_muscle'}},...
                                                'do_plot',false));
td = addCorrelatedNoise(td,struct('signals',{{'linmodel_S1_muscle'}},'noise_covar',muscle_noise_covar));
td = getPCA(td,struct('signals',{{'linmodel_S1_muscle_noisy'}}));

% bootstrap!
n_boot = 1000;
% use actual model for this
bootsep_true = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'S1_pca'}},'mdl',actual_mdl)),td');

bootsep_handle = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_handle_pca',1:4}})),td');
bootsep_handle_actual = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_handle_pca'}},'mdl',actual_mdl)),td');
bootsep_handle_actual_noisy = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_handle_noisy_pca'}},'mdl',actual_mdl)),td');

bootsep_muscle = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_muscle_pca',1:length(opensim_idx)}})),td');
% bootsep_muscle = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_muscle_pca'}})),td');
bootsep_muscle_actual = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_muscle_pca'}},'mdl',actual_mdl)),td');
bootsep_muscle_actual_noisy = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_muscle_noisy_pca'}},'mdl',actual_mdl)),td');

% get separability CIs
sep_true = prctile(bootsep_true,[2.5 97.5]);

sep_handle = prctile(bootsep_handle,[2.5 97.5]);
sep_handle_actual = prctile(bootsep_handle_actual,[2.5 97.5]);
sep_handle_actual_noisy = prctile(bootsep_handle_actual_noisy,[2.5 97.5]);

sep_muscle = prctile(bootsep_muscle,[2.5 97.5]);
sep_muscle_actual = prctile(bootsep_muscle_actual,[2.5 97.5]);
sep_muscle_actual_noisy = prctile(bootsep_muscle_actual_noisy,[2.5 97.5]);

% plot separabilities with CIs
figure
barh(fliplr([mean(bootsep_true) mean(bootsep_handle) mean(bootsep_handle_actual) mean(bootsep_handle_actual_noisy)...
    mean(bootsep_muscle) mean(bootsep_muscle_actual) mean(bootsep_muscle_actual_noisy)]))
shading flat
hold on
plot(fliplr([sep_true' sep_handle' sep_handle_actual' sep_handle_actual_noisy'...
                sep_muscle' sep_muscle_actual' sep_muscle_actual_noisy']),[1:7;1:7],'-k','linewidth',2)
plot([0.5 0.5],[0 8],'--k','linewidth',2)
set(gca,'box','off','tickdir','out','yticklabel',fliplr({'Actual','Handle Vel/Force','Handle Vel/Force w/Actual Disc',...
                                                    'Noisy Handle Vel/Force w/Actual Disc','Muscle Vel',...
                                                    'Muscle Vel w/Actual Disc','Noisy Muscle Vel w/Actual Disc'}),...
                                                    'xtick',[0 0.5 1])

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
