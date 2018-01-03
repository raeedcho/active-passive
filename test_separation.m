%% Get PCA for act vs pas
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
% clean nans out...?
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
td_act = binTD(td_act,15);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_pas = binTD(td_pas,15);
td = cat(2,td_act,td_pas);

%% set figure savename
figsave_name = [td(1).monkey '_' datestr(datenum(td(1).date),'yyyymmdd') '_'];

%% get splits on array
if exist('ant_chans','var')
    % array_signals = find(ismember(td(1).S1_unit_guide(:,1),ant_chans));
    array_signals = find(~ismember(td(1).S1_unit_guide(:,1),ant_chans));
    % array_signals = sort(randperm(length(td(1).S1_unit_guide(:,1)),floor(length(td(1).S1_unit_guide(:,1))/2))');
    % array_signals = 1:length(td(1).S1_unit_guide);
else
    array_signals = 1:length(td(1).S1_unit_guide);
end

%% test separability for actual spikes
td = sqrtTransform(td,'S1_spikes');
% td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
td = getPCA(td,struct('signals',{{'S1_spikes',array_signals}}));

figure(1)
[actual_sep,actual_mdl] = test_sep(td,struct('signals',{{'S1_pca'}},'do_plot',true));
title(['Actual separability - ' num2str(actual_sep)])

% then for directional separability/other view
[~,td_act] = getTDidx(td,'ctrHoldBump',false);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
signal_act = cat(1,td_act.S1_pca);
signal_pas = cat(1,td_pas.S1_pca);
% plot active as filled, passive as open
bump_colors = linspecer(4);
act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;
figure(2)
hold all
scatter3(signal_act(:,1),signal_act(:,2),signal_act(:,3),50,bump_colors(act_dir_idx,:),'filled')
scatter3(signal_pas(:,1),signal_pas(:,2),signal_pas(:,3),100,bump_colors(pas_dir_idx,:),'o','linewidth',2)
ylim = get(gca,'ylim');
zlim = get(gca,'zlim');
set(gca,'box','off','tickdir','out')
axis equal
axis off
title 'directional sep'

%% Try fabricating trial_data with linear models based on handle kinematics and force
% get models for force and velocity from actpas data
[td,model_info] = getModel(td,struct('model_type','linmodel',...
    'model_name','S1_handle','in_signals',{{'vel';'force'}},...
    'out_signals',{{'S1_spikes',array_signals}}));

td = getPCA(td,struct('signals',{{'linmodel_S1_handle'}}));
% get correlated noise model
handle_noise_covar = getNoiseCovar(td,struct('actual_signals',{{'S1_spikes',array_signals}},'modeled_signals',{{'linmodel_S1_handle'}},...
                                                'do_plot',false));
td = addCorrelatedNoise(td,struct('signals',{{'linmodel_S1_handle'}},'noise_covar',handle_noise_covar));
td = getPCA(td,struct('signals',{{'linmodel_S1_handle_noisy'}}));

figure(3)
[velforce_sep,velforce_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_handle_pca',1:4}},'do_plot',true));
title(['Velforce separability - ' num2str(velforce_sep)])
figure(4)
[velforce_sep,velforce_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_handle_pca'}},'mdl',actual_mdl,'do_plot',true));
title(['Velforce w/ actual disc. separability - ' num2str(velforce_sep)])
figure(5)
[velforce_sep,velforce_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_handle_noisy_pca'}},'mdl',actual_mdl,'do_plot',true));
title(['Noisy velforce w/ actual disc. separability - ' num2str(velforce_sep)])

% % then for directional separability/other view
% [~,td_act] = getTDidx(td,'ctrHoldBump',false);
% [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
% signal_act = cat(1,td_act.linmodel_S1_handle_pca);
% signal_pas = cat(1,td_pas.linmodel_S1_handle_pca);
% % plot active as filled, passive as open
% bump_colors = linspecer(4);
% act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
% pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;
% subplot(2,3,2)
% hold all
% scatter3(signal_act(:,1),signal_act(:,2),signal_act(:,3),50,bump_colors(act_dir_idx,:),'filled')
% scatter3(signal_pas(:,1),signal_pas(:,2),signal_pas(:,3),100,bump_colors(pas_dir_idx,:),'o','linewidth',2)
% ylim = get(gca,'ylim');
% zlim = get(gca,'zlim');
% set(gca,'box','off','tickdir','out')
% axis equal
% axis off

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
    'out_signals',{{'S1_spikes',array_signals}}));

td = getPCA(td,struct('signals',{{'linmodel_S1_muscle'}}));
% get correlated noise model
muscle_noise_covar = getNoiseCovar(td,struct('actual_signals',{{'S1_spikes',array_signals}},'modeled_signals',{{'linmodel_S1_muscle'}},...
                                                'do_plot',false));
td = addCorrelatedNoise(td,struct('signals',{{'linmodel_S1_muscle'}},'noise_covar',muscle_noise_covar));
td = getPCA(td,struct('signals',{{'linmodel_S1_muscle_noisy'}}));

figure(6)
% [muscle_sep,muscle_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_muscle_pca',1:length(opensim_idx)}},'do_plot',true));
[muscle_sep,muscle_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_muscle_pca'}},'do_plot',true));
title(['Muscle separability - ' num2str(muscle_sep)])
figure(7)
[muscle_sep,muscle_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_muscle_pca'}},'mdl',actual_mdl,'do_plot',true));
title(['Muscle w/ actual disc. separability - ' num2str(muscle_sep)])
figure(8)
[muscle_sep,muscle_mdl] = test_sep(td,struct('signals',{{'linmodel_S1_muscle_noisy_pca'}},'mdl',actual_mdl,'do_plot',true));
title(['Noisy muscle w/ actual disc. separability - ' num2str(muscle_sep)])

% % then for directional separability/other view
% [~,td_act] = getTDidx(td,'ctrHoldBump',false);
% [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
% signal_act = cat(1,td_act.linmodel_S1_muscle_pca);
% signal_pas = cat(1,td_pas.linmodel_S1_muscle_pca);
% % plot active as filled, passive as open
% bump_colors = linspecer(4);
% act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
% pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;
% subplot(2,3,3)
% hold all
% scatter3(signal_act(:,1),signal_act(:,2),signal_act(:,3),50,bump_colors(act_dir_idx,:),'filled')
% scatter3(signal_pas(:,1),signal_pas(:,2),signal_pas(:,3),100,bump_colors(pas_dir_idx,:),'o','linewidth',2)
% ylim = get(gca,'ylim');
% zlim = get(gca,'zlim');
% set(gca,'box','off','tickdir','out')
% axis equal
% axis off

%% get boostrapped separability values

% get correlated noise models

% bootstrap!
n_boot = 1000;
% use actual model for this
bootsep_true = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'S1_pca'}},'mdl',actual_mdl)),td');

bootsep_handle = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_handle_pca',1:4}})),td');
bootsep_handle_actual = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_handle_pca'}},'mdl',actual_mdl)),td');
bootsep_actual_handle = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_handle_noisy_pca'}},'mdl',velforce_mdl)),td');

% bootsep_muscle = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_muscle_pca',1:length(opensim_idx)}})),td');
bootsep_muscle = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_muscle_pca'}})),td');
bootsep_muscle_actual = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_muscle_pca'}},'mdl',actual_mdl)),td');
bootsep_actual_muscle = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'linmodel_S1_muscle_noisy_pca'}},'mdl',muscle_mdl)),td');

% get separability CIs
sep_true = prctile(bootsep_true,[2.5 97.5]);

sep_handle = prctile(bootsep_handle,[2.5 97.5]);
sep_handle_actual = prctile(bootsep_handle_actual,[2.5 97.5]);
sep_actual_handle = prctile(bootsep_actual_handle,[2.5 97.5]);

sep_muscle = prctile(bootsep_muscle,[2.5 97.5]);
sep_muscle_actual = prctile(bootsep_muscle_actual,[2.5 97.5]);
sep_actual_muscle = prctile(bootsep_actual_muscle,[2.5 97.5]);

% plot separabilities with CIs
figure
barh(fliplr([mean(bootsep_true) mean(bootsep_handle) mean(bootsep_handle_actual) mean(bootsep_actual_handle)...
    mean(bootsep_muscle) mean(bootsep_muscle_actual) mean(bootsep_actual_muscle)]))
shading flat
hold on
plot(fliplr([sep_true' sep_handle' sep_handle_actual' sep_actual_handle'...
                sep_muscle' sep_muscle_actual' sep_actual_muscle']),[1:7;1:7],'-k','linewidth',2)
plot([0.5 0.5],[0 8],'--k','linewidth',2)
set(gca,'box','off','tickdir','out','yticklabel',fliplr({'Actual','Handle Vel/Force','Handle Vel/Force w/Actual Disc',...
                                                    'Actual w/Handle Vel/Force Disc','Muscle Vel',...
                                                    'Muscle Vel w/Actual Disc','Actual w/Muscle Vel Disc'}),...
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
