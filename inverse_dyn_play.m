%% Extract trial data from cds files
params.array_alias = {'LeftS1Area2','S1'};
% params.exclude_units = [255];
params.event_list = {'ctrHoldBump';'bumpTime';'bumpDir'};
params.trial_results = {'R','A','F','I'};
meta = struct('task','COactpas');
params.meta = meta;
trial_data_actpas = parseFileByTrial(cds,params);

%% plot some inverse dynamics stuff
[~,td] = getTDidx(trial_data_actpas,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false,'target_direction',pi/2);
td_act = trimTD(td_act,{'idx_movement_on',-10},{'idx_movement_on',50});

[~,td_pas] = getTDidx(td,'ctrHoldBump',true,'bumpDir',90);
td_pas = trimTD(td_pas,{'idx_bumpTime',-10},{'idx_bumpTime',50});

% td_still = td;
% 
% % check if any events happen within 500 ms of start time
% 
% % otherwise just average
% td_still = trimTD(td_still,{'idx_startTime',0},{'idx_startTime',50});
% td_still = trialAverage(td_still,{'ctrHoldBump'});

% clean nans out...?
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);

% Trial average
% td_act = trialAverage(td_act,{'target_direction'});
% td_pas = trialAverage(td_pas,{'bumpDir'});

% plot active as filled, passive as open
% bump_colors = linspecer(4);

% plot dynamics of shoulder flexion
% trial_num = 3;


for trial_num = 1:30

    flex_vel_idx = find(strcmp('shoulder_flexion_vel',td(trial_num).opensim_names));
    flex_mom_idx = find(strcmp('shoulder_flexion_moment',td(trial_num).opensim_names));

    % plot angular velocity and moment
    figure(1234)
    clf;
    subplot 231
    yyaxis left

    times_act = -10:length(td_act(trial_num).opensim(:,flex_vel_idx))-11;
    times_act = times_act*td_act(trial_num).bin_size;
    %calculate velocity
    %         joint_acc = gradient(td_act(trial_num).opensim(:,3),1);
    joint_vel = td_act(trial_num).opensim(:,flex_vel_idx);
    joint_mom = td_act(trial_num).opensim(:,flex_mom_idx);
    plot(times_act,joint_vel,'linewidth',2)
    yyaxis right
    plot(times_act,joint_mom,'linewidth',2)
    hold on
    ylims = get(gca,'ylim');
    plot(times_act(repmat(td_act(trial_num).idx_movement_on,2,1)),ylims,'k--','linewidth',2)
    title 'Active shoulder flexion'
    set(gca,'box','off','tickdir','out')

    subplot 234
    yyaxis left

    times_pas = -10:length(td_pas(trial_num).opensim(:,flex_vel_idx))-11;
    times_pas = times_pas*td_pas(trial_num).bin_size;
    %calculate velocity
    joint_vel = td_pas(trial_num).opensim(:,flex_vel_idx);
    joint_mom = td_pas(trial_num).opensim(:,flex_mom_idx);
    plot(times_pas,joint_vel,'linewidth',2)
    yyaxis right
    plot(times_pas,joint_mom,'linewidth',2)
    hold on
    ylims = get(gca,'ylim');
    plot(times_pas(repmat(td_pas(trial_num).idx_bumpTime,2,1)),ylims,'k--','linewidth',2)
    title 'Passive shoulder flexion'
    set(gca,'box','off','tickdir','out')

    xlabel 'Time post-bump/move (s)'

    % plot power
    subplot 232

    times_act = -10:length(td_act(trial_num).opensim(:,flex_vel_idx))-11;
    times_act = times_act*td_act(trial_num).bin_size;
    %calculate velocity
    %         joint_acc = gradient(td_act(trial_num).opensim(:,3),1);
    joint_vel = td_act(trial_num).opensim(:,flex_vel_idx);
    joint_mom = td_act(trial_num).opensim(:,flex_mom_idx);
    plot(times_act,joint_vel.*joint_mom,'linewidth',2)
    hold on
    ylims = get(gca,'ylim');
    plot(times_act(repmat(td_act(trial_num).idx_movement_on,2,1)),ylims,'k--','linewidth',2)
    title 'Active shoulder power'
    set(gca,'box','off','tickdir','out')

    subplot 235

    times_pas = -10:length(td_pas(trial_num).opensim(:,flex_vel_idx))-11;
    times_pas = times_pas*td_pas(trial_num).bin_size;
    %calculate velocity
    joint_vel = td_pas(trial_num).opensim(:,flex_vel_idx);
    joint_mom = td_pas(trial_num).opensim(:,flex_mom_idx);
    plot(times_pas,joint_vel.*joint_mom,'linewidth',2)
    hold on
    ylims = get(gca,'ylim');
    plot(times_pas(repmat(td_pas(trial_num).idx_bumpTime,2,1)),ylims,'k--','linewidth',2)
    title 'Passive shoulder power'
    set(gca,'box','off','tickdir','out')

    xlabel 'Time post-bump/move (s)'

    % plot handle power
    subplot 233

    times_act = -10:length(td_act(trial_num).opensim(:,flex_vel_idx))-11;
    times_act = times_act*td_act(trial_num).bin_size;
    %calculate velocity
    %         joint_acc = gradient(td_act(trial_num).opensim(:,3),1);
    joint_vel = td_act(trial_num).vel;
    joint_mom = td_act(trial_num).force(:,1:2);
    plot(times_act,sum(joint_vel.*joint_mom,2),'linewidth',2)
    hold on
    ylims = get(gca,'ylim');
    plot(times_act(repmat(td_act(trial_num).idx_movement_on,2,1)),ylims,'k--','linewidth',2)
    title 'Active handle power'
    set(gca,'box','off','tickdir','out')

    subplot 236

    times_pas = -10:length(td_pas(trial_num).opensim(:,flex_vel_idx))-11;
    times_pas = times_pas*td_pas(trial_num).bin_size;
    %calculate velocity
    joint_vel = td_pas(trial_num).vel;
    joint_mom = td_pas(trial_num).force(:,1:2);
    plot(times_pas,sum(joint_vel.*joint_mom,2),'linewidth',2)
    hold on
    ylims = get(gca,'ylim');
    plot(times_pas(repmat(td_pas(trial_num).idx_bumpTime,2,1)),ylims,'k--','linewidth',2)
    title 'Passive handle power'
    set(gca,'box','off','tickdir','out')

    xlabel 'Time post-bump/move (s)'

    waitforbuttonpress;
end

