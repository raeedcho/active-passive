%% Look at forces and speeds throughout movement
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false,'target_direction',pi/2);
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});

[~,td_pas] = getTDidx(td,'ctrHoldBump',true,'bumpDir',90);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});

% clean nans out...?
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);

% Trial average
td_act = trialAverage(td_act,{'target_direction'});
td_pas = trialAverage(td_pas,{'bumpDir'});

figure
for trial_num = 1:length(td_act)
    %calculate acceleration
    act_speed = sqrt(sum(td_act(trial_num).vel.^2,2));
    
    plot(act_speed,'-b','linewidth',1)
    hold on
    
    % ylims = get(gca,'ylim');
    % plot(repmat(td_act(trial_num).idx_movement_on,2,1),ylims,'r--','linewidth',1)
    % plot(repmat(td_act(trial_num).idx_peak_speed,2,1),ylims,'g--','linewidth',1)
    
    set(gca,'box','off','tickdir','out')
end

for trial_num = 1:length(td_pas)
    %calculate acceleration
    pas_speed = sqrt(sum(td_pas(trial_num).vel.^2,2));
    
    plot(pas_speed,'-g','linewidth',1)
    hold on
    
    % ylims = get(gca,'ylim');
    % plot(repmat(td_pas(trial_num).idx_bumpTime,2,1),ylims,'k--','linewidth',1)
    
    set(gca,'box','off','tickdir','out')
end
