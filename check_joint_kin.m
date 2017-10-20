%% Get PCA for act vs pas
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);

% get still handle data (no word for start of center hold)
minCH = min(cat(1,td_act.ctrHold));
bin_size = td_act(1).bin_size;
still_bins = floor(minCH/bin_size);
td_still = trimTD(td_act,{'idx_goCueTime',-still_bins},{'idx_goCueTime',0});

% Get td_act and td_pas
num_bins_before = floor(still_bins/2);
num_bins_after = 30;

td_act = trimTD(td_act,{'idx_movement_on',-num_bins_before},{'idx_movement_on',num_bins_after});
td_pas = trimTD(td_pas,{'idx_bumpTime',-num_bins_before},{'idx_bumpTime',num_bins_after});

[td_act_avg,act_cond_idx] = trialAverage(td_act,'target_direction');
[td_pas_avg,pas_cond_idx] = trialAverage(td_pas,'bumpDir');

times = -num_bins_before:num_bins_after;
times = times*td_act_avg(1).bin_size;
opensim_idx = find(contains(td(1).opensim_names,'_vel'));
subplot_idx = [6 2 4 8];
colors = linspecer(4);
for i = 1:numel(opensim_idx)
    figure

    % plot to figure out limits for uniform scale
    for j = 1:4
        plot(times,td_act_avg(j).opensim(:,opensim_idx(i)),'-','linewidth',2,'color',colors(j,:))
        hold on
        plot(times,td_pas_avg(j).opensim(:,opensim_idx(i)),'--','linewidth',2,'color',colors(j,:))
        set(gca,'box','off','tickdir','out')
    end

    % now plot for real
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    clf
    for j = 1:4
        subplot(3,3,subplot_idx(j))
        plot(times,td_act_avg(j).opensim(:,opensim_idx(i)),'-','linewidth',2,'color',colors(j,:))
        hold on
        plot(times,td_pas_avg(j).opensim(:,opensim_idx(i)),'--','linewidth',2,'color',colors(j,:))

        % get CIs
        temp_td = td_act(act_cond_idx{j});
        opensim_data = cat(3,temp_td.opensim);
        opensimCI = prctile(opensim_data,[10 90],3);
        opensimCIlow = squeeze(opensimCI(:,:,1));
        opensimCIhigh = squeeze(opensimCI(:,:,2));
        opensim_patch = cat(1,opensimCIlow,flipud(opensimCIhigh));
        times_patch = cat(2,times,fliplr(times));
        h = patch(times_patch',opensim_patch(:,opensim_idx(i)),colors(j,:));
        set(h,'edgecolor','none','facealpha',0.5)

        temp_td = td_pas(pas_cond_idx{j});
        opensim_data = cat(3,temp_td.opensim);
        opensimCI = prctile(opensim_data,[10 90],3);
        opensimCIlow = squeeze(opensimCI(:,:,1));
        opensimCIhigh = squeeze(opensimCI(:,:,2));
        opensim_patch = cat(1,opensimCIlow,flipud(opensimCIhigh));
        times_patch = cat(2,times,fliplr(times));
        h = patch(times_patch',opensim_patch(:,opensim_idx(i)),colors(j,:));
        set(h,'edgecolor','none','facealpha',0.5)


        plot([0 0],ylim,'--k','linewidth',3)
        set(gca,'box','off','tickdir','out','xlim',xlim,'ylim',ylim)
    end
    subplot(3,3,2)
    title(['Angular velocity of ' td(1).opensim_names(opensim_idx(i)) ' (deg/s)'])
end

% now do straight handle velocity
vel_names = {'x','y'};
for i = 1:2
    figure

    % plot to figure out limits for uniform scale
    for j = 1:4
        plot(times,td_act_avg(j).vel(:,i),'-','linewidth',2,'color',colors(j,:))
        hold on
        plot(times,td_pas_avg(j).vel(:,i),'--','linewidth',2,'color',colors(j,:))
        set(gca,'box','off','tickdir','out')
    end

    % now plot for real
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    clf
    for j = 1:4
        subplot(3,3,subplot_idx(j))
        plot(times,td_act_avg(j).vel(:,i),'-','linewidth',2,'color',colors(j,:))
        hold on
        plot(times,td_pas_avg(j).vel(:,i),'--','linewidth',2,'color',colors(j,:))
        plot([0 0],ylim,'--k','linewidth',3)
        set(gca,'box','off','tickdir','out','xlim',xlim,'ylim',ylim)

        % get CIs
        temp_td = td_act(act_cond_idx{j});
        opensim_data = cat(3,temp_td.vel);
        opensimCI = prctile(opensim_data,[10 90],3);
        opensimCIlow = squeeze(opensimCI(:,:,1));
        opensimCIhigh = squeeze(opensimCI(:,:,2));
        opensim_patch = cat(1,opensimCIlow,flipud(opensimCIhigh));
        times_patch = cat(2,times,fliplr(times));
        h = patch(times_patch',opensim_patch(:,i),colors(j,:));
        set(h,'edgecolor','none','facealpha',0.5)

        temp_td = td_pas(pas_cond_idx{j});
        opensim_data = cat(3,temp_td.vel);
        opensimCI = prctile(opensim_data,[10 90],3);
        opensimCIlow = squeeze(opensimCI(:,:,1));
        opensimCIhigh = squeeze(opensimCI(:,:,2));
        opensim_patch = cat(1,opensimCIlow,flipud(opensimCIhigh));
        times_patch = cat(2,times,fliplr(times));
        h = patch(times_patch',opensim_patch(:,i),colors(j,:));
        set(h,'edgecolor','none','facealpha',0.5)
    end
    subplot(3,3,2)
    title(['Handle velocity of ' vel_names(i) ' (cm/s)'])
end
