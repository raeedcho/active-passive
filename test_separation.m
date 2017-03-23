%%
[~,td] = getTDidx(trial_data_actpas,'result','R');
td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime'));

%% Get PCA for act vs pas

[~,td] = getTDidx(trial_data_actpas,'result','R');
td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime'));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
% td_act = trimTD(td_act,{'idx_goCueTime',35},{'idx_goCueTime',50});
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
td_act = binTD(td_act,15);

[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_pas = binTD(td_pas,15);
td = cat(2,td_act,td_pas);

% td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
td = sqrtTransform(td,'S1_spikes');
[td,pca_info] = getPCA(td,struct('signals',{'S1_spikes'}));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);

% clean nans out...?
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);

% plot active as filled, passive as open
bump_colors = linspecer(4);
S1_pca_act = cat(1,td_act.S1_pca);
S1_pca_pas = cat(1,td_pas.S1_pca);
act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;

% plot dots
figure
hold all
scatter3(S1_pca_act(:,1),S1_pca_act(:,2),S1_pca_act(:,3),50,bump_colors(act_dir_idx,:),'filled')
scatter3(S1_pca_pas(:,1),S1_pca_pas(:,2),S1_pca_pas(:,3),50,bump_colors(pas_dir_idx,:))
axis equal

%% get PCA traces
[~,td] = getTDidx(trial_data_actpas,'result','R');

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_goCueTime',35},{'idx_goCueTime',50});

[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
td = cat(2,td_act,td_pas);

td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
% td = sqrtTransform(td,'S1_spikes');
td = getPCA(td,struct('signals',{'S1_spikes'}));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);

% clean nans out...?
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);

% Trial average
% td_act = trialAverage(td_act,{'target_direction'});
% td_pas = trialAverage(td_pas,{'bumpDir'});

% plot active as filled, passive as open
bump_colors = linspecer(4);

% plot dots
figure
hold all
for trial_ctr = 1:length(td_act)
    S1_pca_act = td_act(trial_ctr).S1_pca;
    act_dir_idx = floor(td_act(trial_ctr).target_direction/(pi/2))+1;
    plot3(S1_pca_act(:,1),  S1_pca_act(:,2),  S1_pca_act(:,3),  '-','Color',bump_colors(act_dir_idx,:),'linewidth',2)
    plot3(S1_pca_act(end,1),S1_pca_act(end,2),S1_pca_act(end,3),'.','Color',bump_colors(act_dir_idx,:),'markersize',30)
end
for trial_ctr = 1:length(td_pas)
    S1_pca_pas = td_pas(trial_ctr).S1_pca;
    pas_dir_idx = floor(td_pas(trial_ctr).bumpDir/(90))+1;
    plot3(S1_pca_pas(:,1),  S1_pca_pas(:,2),  S1_pca_pas(:,3),  '--','Color',bump_colors(pas_dir_idx,:),'linewidth',2)
    plot3(S1_pca_pas(end,1),S1_pca_pas(end,2),S1_pca_pas(end,3),'o', 'Color',bump_colors(pas_dir_idx,:),'markersize',10,'linewidth',2)
end
axis equal



