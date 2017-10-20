
% %% Look at forces and speeds throughout movement
% [~,td] = getTDidx(trial_data_actpas,'result','R');
% 
% td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
% 
% [~,td_act] = getTDidx(td,'ctrHoldBump',false,'target_direction',pi/2);
% td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',12});
% 
% [~,td_pas] = getTDidx(td,'ctrHoldBump',true,'bumpDir',90);
% td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',12});
% 
% % td_still = td;
% % 
% % % check if any events happen within 500 ms of start time
% % 
% % % otherwise just average
% % td_still = trimTD(td_still,{'idx_startTime',0},{'idx_startTime',50});
% % td_still = trialAverage(td_still,{'ctrHoldBump'});
% 
% % clean nans out...?
% nanners = isnan(cat(1,td_act.target_direction));
% td_act = td_act(~nanners);
% 
% Trial average
% td_act = trialAverage(td_act,{'target_direction'});
% td_pas = trialAverage(td_pas,{'bumpDir'});

% plot active as filled, passive as open
% bump_colors = linspecer(4);

% plot dynamics of shoulder flexion
% trial_num = 3;

% figure
% for trial_num = 1:length(td_act)
%     %calculate acceleration
%     act_speed = sqrt(sum(td_act(trial_num).vel.^2,2));
%     act_force = sqrt(sum(td_act(trial_num).force.^2,2));
%     
%     subplot 211
%     
%     yyaxis left
%     plot(act_speed,'-','linewidth',1)
%     hold on
%     
%     yyaxis right
%     plot(act_force,'-','linewidth',1)
%     hold on
%     
%     ylims = get(gca,'ylim');
%     % plot(repmat(td_act(trial_num).idx_goCueTime,2,1),ylims,'k--','linewidth',1)
%     plot(repmat(td_act(trial_num).idx_movement_on,2,1),ylims,'r--','linewidth',1)
%     plot(repmat(td_act(trial_num).idx_peak_speed,2,1),ylims,'g--','linewidth',1)
%     
%     
%     title 'Active movements'
%     set(gca,'box','off','tickdir','out')
% end
% 
% for trial_num = 1:length(td_pas)
%     %calculate acceleration
%     pas_speed = sqrt(sum(td_pas(trial_num).vel.^2,2));
%     pas_force = sqrt(sum(td_pas(trial_num).force.^2,2));
%     
%     subplot 212
%     
%     yyaxis left
%     plot(pas_speed,'-','linewidth',1)
%     hold on
%     
%     yyaxis right
%     plot(pas_force,'-','linewidth',1)
%     hold on
%     
%     ylims = get(gca,'ylim');
%     plot(repmat(td_pas(trial_num).idx_bumpTime,2,1),ylims,'k--','linewidth',1)
%     
%     title 'Passive movements'
%     set(gca,'box','off','tickdir','out')
% end

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

td = sqrtTransform(td,'S1_spikes');
% td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
td = getPCA(td,struct('signals',{{'S1_spikes'}}));
test_sep(td,struct('signals',{{'S1_pca'}},'do_plot',true))

%% Try fabricating trial_data with linear models based on behavior
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
td_act = binTD(td_act,15);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_pas = binTD(td_pas,15);
td = cat(2,td_act,td_pas);

% get models for force and velocity from actpas data
[td,model_info] = getModel(td,struct('model_type','linmodel',...
    'model_name','S1','in_signals',{{'vel';'force'}},...
    'out_signals',{'S1_spikes'}));

td = getPCA(td,struct('signals',{{'linmodel_S1'}}));
test_sep(td,struct('signals',{{'linmodel_S1_pca',1:4}},'do_plot',true))

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

%% Try fabricating trial_data with linear models based on joint stuff
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
td_act = binTD(td_act,15);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_pas = binTD(td_pas,15);
td = cat(2,td_act,td_pas);

% get models for force and velocity from actpas data
% opensim_idx = find(contains(td(1).opensim_names,'_moment'));
opensim_idx = find(contains(td(1).opensim_names,'_vel'));
% opensim_idx = find(contains(td(1).opensim_names,'_vel') & ~contains(td(1).opensim_names,'wrist') & ~contains(td(1).opensim_names,'radial'));
% opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
% opensim_idx = find(contains(td(1).opensim_names,'_vel') | contains(td(1).opensim_names,'_moment'));
% opensim_idx = find( (contains(td(1).opensim_names,'_vel') | contains(td(1).opensim_names,'_moment')) & ~contains(td(1).opensim_names,'wrist') & ~contains(td(1).opensim_names,'radial'));
[td,model_info] = getModel(td,struct('model_type','linmodel',...
    'model_name','S1','in_signals',...
    {{'opensim',opensim_idx}},...
    'out_signals',{'S1_spikes'}));

td = getPCA(td,struct('signals',{{'linmodel_S1'}}));
test_sep(td,struct('signals',{{'linmodel_S1_pca',1:length(opensim_idx)}},'do_plot',true))

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
