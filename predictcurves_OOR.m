%% Use intended direction and force
[~,td] = getTDidx(trial_data,'result','R');
td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
td = trimTD(td,{'idx_movement_on',0},{'idx_endTargHoldTime',0});

% fix 360 force direction
for i = 1:length(td)
    td(i).force_direction = mod(td(i).force_direction,360);
end

% get actual firing rates
for i = 1:length(td)
    td(i).S1_spikes = td(i).S1_spikes/td(i).bin_size;
end

% create models
[train_idx,test_idx] = crossvalind('HoldOut',length(td),0.6);
[td,model_info] = getModel(td,struct('model_type','linmodel',...
    'model_name','S1_handle','in_signals',{{'vel'}},...
    'out_signals',{'S1_spikes'},'train_idx',find(train_idx)));
[td,model_info] = getModel(td,struct('model_type','linmodel',...
    'model_name','S1_velforce','in_signals',{{'vel';'force'}},...
    'out_signals',{'S1_spikes'},'train_idx',find(train_idx)));
opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
[td,model_info] = getModel(td,struct('model_type','linmodel',...
    'model_name','S1_muscle','in_signals', {{'opensim',opensim_idx}},...
    'out_signals',{'S1_spikes'},'train_idx',find(train_idx)));

% get PDs for each model
n_boot = 1000;
pd_true = get_targ_PDs(td,struct('out_signals','S1_spikes','out_signal_names',{td(1).S1_unit_guide},'num_boots',n_boot,'trial_idx',find(test_idx)));
pd_velforce = get_targ_PDs(td,struct('out_signals','linmodel_S1_velforce','out_signal_names',{td(1).S1_unit_guide},'num_boots',n_boot,'trial_idx',find(test_idx)));
pd_handle = get_targ_PDs(td,struct('out_signals','linmodel_S1_handle','out_signal_names',{td(1).S1_unit_guide},'num_boots',n_boot,'trial_idx',find(test_idx)));
pd_muscle = get_targ_PDs(td,struct('out_signals','linmodel_S1_muscle','out_signal_names',{td(1).S1_unit_guide},'num_boots',n_boot,'trial_idx',find(test_idx)));

% create confusion matrix for tuning
CIthresh = pi/2;
isTuned_targ_true = checkIsTuned(pd_true,struct('move_corr','target_direction','CIthresh',CIthresh));
isTuned_load_true = checkIsTuned(pd_true,struct('move_corr','force_direction','CIthresh',CIthresh));
isTuned_targ_muscle = checkIsTuned(pd_muscle,struct('move_corr','target_direction','CIthresh',CIthresh));
isTuned_load_muscle = checkIsTuned(pd_muscle,struct('move_corr','force_direction','CIthresh',CIthresh));
isTuned_targ_handle = checkIsTuned(pd_handle,struct('move_corr','target_direction','CIthresh',CIthresh));
isTuned_load_handle = checkIsTuned(pd_handle,struct('move_corr','force_direction','CIthresh',CIthresh));
isTuned_targ_velforce = checkIsTuned(pd_velforce,struct('move_corr','target_direction','CIthresh',CIthresh));
isTuned_load_velforce = checkIsTuned(pd_velforce,struct('move_corr','force_direction','CIthresh',CIthresh));
num_units = size(td(1).S1_unit_guide,1);

targ_confusion = [sum(isTuned_targ_true & isTuned_targ_muscle) sum(isTuned_targ_true & ~isTuned_targ_muscle);...
                 sum(~isTuned_targ_true & isTuned_targ_muscle) sum(~isTuned_targ_true & ~isTuned_targ_muscle)]/num_units;
load_confusion = [sum(isTuned_load_true & isTuned_load_muscle) sum(isTuned_load_true & ~isTuned_load_muscle);...
                 sum(~isTuned_load_true & isTuned_load_muscle) sum(~isTuned_load_true & ~isTuned_load_muscle)]/num_units;
targ_confusion_velforce = [sum(isTuned_targ_true & isTuned_targ_velforce) sum(isTuned_targ_true & ~isTuned_targ_velforce);...
                 sum(~isTuned_targ_true & isTuned_targ_velforce) sum(~isTuned_targ_true & ~isTuned_targ_velforce)]/num_units;
load_confusion_velforce = [sum(isTuned_load_true & isTuned_load_velforce) sum(isTuned_load_true & ~isTuned_load_velforce);...
                 sum(~isTuned_load_true & isTuned_load_velforce) sum(~isTuned_load_true & ~isTuned_load_velforce)]/num_units;
targ_confusion_handle = [sum(isTuned_targ_true & isTuned_targ_handle) sum(isTuned_targ_true & ~isTuned_targ_handle);...
                 sum(~isTuned_targ_true & isTuned_targ_handle) sum(~isTuned_targ_true & ~isTuned_targ_handle)]/num_units;
load_confusion_handle = [sum(isTuned_load_true & isTuned_load_handle) sum(isTuned_load_true & ~isTuned_load_handle);...
                 sum(~isTuned_load_true & isTuned_load_handle) sum(~isTuned_load_true & ~isTuned_load_handle)]/num_units;

%% plots?
% plot muscle active tuning against passive tuning
h1 = figure;
subplot(1,2,2)
isTuned_load = isTuned_load_true & isTuned_load_velforce;
comparePDs(pd_true(isTuned_load,:),pd_velforce(isTuned_load,:),struct('move_corr','force_direction'),'ro','linewidth',2)
hold on
isTuned_load = isTuned_load_true & isTuned_load_muscle;
comparePDs(pd_true(isTuned_load,:),pd_muscle(isTuned_load,:),struct('move_corr','force_direction'),'co','linewidth',2)
% isTuned_load = isTuned_load_true & isTuned_load_handle;
% comparePDs(pd_true(isTuned_load,:),pd_handle(isTuned_load,:),struct('move_corr','force_direction'),'go','linewidth',2)
plot([-pi pi],[-pi pi],'--k','linewidth',2)
xlabel('Actual load PD')
ylabel('Predicted load PD')
title('Load PDs')
set(get(gca,'ylabel'),'rotation',0,'horizontalalignment','right')
set(gca,'box','off','tickdir','out','xlim',[-pi pi],'ylim',[-pi pi],'xtick',[-pi 0 pi],'ytick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'},'yticklabel',{'-\pi','0','\pi'})
axis equal

subplot(1,2,1)
isTuned_targ = isTuned_targ_true & isTuned_targ_velforce;
comparePDs(pd_true(isTuned_targ,:),pd_velforce(isTuned_targ,:),struct('move_corr','target_direction'),'ro','linewidth',2)
hold on
isTuned_targ = isTuned_targ_true & isTuned_targ_muscle;
comparePDs(pd_true(isTuned_targ,:),pd_muscle(isTuned_targ,:),struct('move_corr','target_direction'),'co','linewidth',2)
% isTuned_targ = isTuned_targ_true & isTuned_targ_handle;
% comparePDs(pd_true(isTuned_targ,:),pd_handle(isTuned_targ,:),struct('move_corr','target_direction'),'go','linewidth',2)
plot([-pi pi],[-pi pi],'--k','linewidth',2)
xlabel('Actual targ PD')
ylabel('Predicted targ PD')
title('Target PDs')
set(get(gca,'ylabel'),'rotation',0,'horizontalalignment','right')
set(gca,'box','off','tickdir','out','xlim',[-pi pi],'ylim',[-pi pi],'xtick',[-pi 0 pi],'ytick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'},'yticklabel',{'-\pi','0','\pi'})
axis equal
