%% plot some inverse dynamics stuff
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);

% get still handle data (no word for start of center hold)
minCH = min(cat(1,td_act.ctrHold));
bin_size = td_act(1).bin_size;
still_bins = floor(minCH/bin_size);
td_still = trimTD(td_act,{'idx_goCueTime',-still_bins},{'idx_goCueTime',0});

% setup
num_bins_before = floor(still_bins/2);
num_bins_after = 30;

% get active movements
td_act = trimTD(td_act,{'idx_movement_on',-num_bins_before},{'idx_movement_on',num_bins_after});

% get passive movements
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',-num_bins_before},{'idx_bumpTime',num_bins_after});

% Trial average
td_act = trialAverage(td_act,{'target_direction'});
td_pas = trialAverage(td_pas,{'bumpDir'});
td_still = trialAverage(td_still,{'ctrHoldBump'});

% set up which index to look at
joint_name = 'shoulder_flexion';
opensim_vel_idx = find(strcmp([joint_name '_vel'],td(1).opensim_names));
opensim_mom_idx = find(strcmp([joint_name '_moment'],td(1).opensim_names));
mom_offset =  mean(td_still.opensim(:,opensim_mom_idx));
% mom_offset = 0;

% for all four directions
for trial_num = 1:4
    % plot angular velocity and moment
    % figure(1234)
    figure
    clf;
    subplot 231
    yyaxis left

    times_act = -num_bins_before:length(td_act(trial_num).opensim(:,opensim_vel_idx))-(num_bins_before+1);
    times_act = times_act*td_act(trial_num).bin_size;
    %calculate velocity
    joint_vel = td_act(trial_num).opensim(:,opensim_vel_idx);
    joint_mom = td_act(trial_num).opensim(:,opensim_mom_idx) - mom_offset;
    plot(times_act,joint_vel,'linewidth',2)
    hold on
    % joint_vel = td_still(trial_num).opensim(:,opensim_vel_idx);
    % plot([times_act(1);times_act(end)],repmat(mean(joint_vel),2,1),'--','linewidth',2)
    yyaxis right
    plot(times_act,joint_mom,'linewidth',2)
    hold on
    % joint_mom = td_still(trial_num).opensim(:,opensim_mom_idx);
    % plot([times_act(1);times_act(end)],repmat(mean(joint_mom),2,1),'--','linewidth',2)
    ylims = get(gca,'ylim');
    plot(times_act(repmat(td_act(trial_num).idx_movement_on,2,1)),ylims,'k--','linewidth',2)
    title(['Active ' joint_name])
    set(gca,'box','off','tickdir','out')

    subplot 234
    yyaxis left

    times_pas = -num_bins_before:length(td_pas(trial_num).opensim(:,opensim_vel_idx))-(num_bins_before+1);
    times_pas = times_pas*td_pas(trial_num).bin_size;
    %calculate velocity
    joint_vel = td_pas(trial_num).opensim(:,opensim_vel_idx);
    joint_mom = td_pas(trial_num).opensim(:,opensim_mom_idx) - mom_offset;
    plot(times_pas,joint_vel,'linewidth',2)
    hold on
    % joint_vel = td_still(trial_num).opensim(:,opensim_vel_idx);
    % plot([times_pas(1);times_pas(end)],repmat(mean(joint_vel),2,1),'--','linewidth',2)
    yyaxis right
    plot(times_pas,joint_mom,'linewidth',2)
    hold on
    % joint_mom = td_still(trial_num).opensim(:,opensim_mom_idx);
    % plot([times_pas(1);times_pas(end)],repmat(mean(joint_mom),2,1),'--','linewidth',2)
    ylims = get(gca,'ylim');
    plot(times_pas(repmat(td_pas(trial_num).idx_bumpTime,2,1)),ylims,'k--','linewidth',2)
    title(['Passive ' joint_name])
    set(gca,'box','off','tickdir','out')

    xlabel 'Time post-bump/move (s)'

    % plot power
    subplot 232

    times_act = -num_bins_before:length(td_act(trial_num).opensim(:,opensim_vel_idx))-(num_bins_before+1);
    times_act = times_act*td_act(trial_num).bin_size;
    %calculate velocity
    %         joint_acc = gradient(td_act(trial_num).opensim(:,3),1);
    joint_vel = td_act(trial_num).opensim(:,opensim_vel_idx);
    joint_mom = td_act(trial_num).opensim(:,opensim_mom_idx) - mom_offset;
    plot(times_act,joint_vel.*joint_mom,'linewidth',2)
    hold on
    ylims = get(gca,'ylim');
    plot(times_act(repmat(td_act(trial_num).idx_movement_on,2,1)),ylims,'k--','linewidth',2)
    title(['Active ' joint_name ' power'])
    set(gca,'box','off','tickdir','out')

    subplot 235

    times_pas = -num_bins_before:length(td_pas(trial_num).opensim(:,opensim_vel_idx))-(num_bins_before+1);
    times_pas = times_pas*td_pas(trial_num).bin_size;
    %calculate velocity
    joint_vel = td_pas(trial_num).opensim(:,opensim_vel_idx);
    joint_mom = td_pas(trial_num).opensim(:,opensim_mom_idx) - mom_offset;
    plot(times_pas,joint_vel.*joint_mom,'linewidth',2)
    hold on
    ylims = get(gca,'ylim');
    plot(times_pas(repmat(td_pas(trial_num).idx_bumpTime,2,1)),ylims,'k--','linewidth',2)
    title(['Passive ' joint_name ' power'])
    set(gca,'box','off','tickdir','out')

    xlabel 'Time post-bump/move (s)'

    % plot handle power
    subplot 233

    times_act = -num_bins_before:length(td_act(trial_num).opensim(:,opensim_vel_idx))-(num_bins_before+1);
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

    times_pas = -num_bins_before:length(td_pas(trial_num).opensim(:,opensim_vel_idx))-(num_bins_before+1);
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

    % waitforbuttonpress;
end
