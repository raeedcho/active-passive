%% Extract trial data from cds files
params.array_alias = {'LeftS1Area2','S1'};
% params.exclude_units = [255];
params.event_list = {'ctrHoldBump';'bumpTime';'bumpDir'};
params.trial_results = {'R','A','F','I'};
meta = struct('task','COactpas');
params.meta = meta;
trial_data_actpas = parseFileByTrial(actpas_cds_speeds,params);

%% plot some inverse dynamics stuff
[~,td] = getTDidx(trial_data_actpas,'result','R');

[~,td_act] = getTDidx(td,'ctrHoldBump',false,'target_direction',pi/2);
td_act = trimTD(td_act,{'idx_startTime',0},{'idx_endTime',0});

[~,td_pas] = getTDidx(td,'ctrHoldBump',true,'bumpDir',90);
td_pas = trimTD(td_pas,{'idx_startTime',0},{'idx_bumpTime',50});

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


for trial_num = 1:15
    figure
    subplot 211
    yyaxis left
    
    %calculate acceleration
    joint_acc = gradient(td_act(trial_num).opensim(:,3),1);
    plot(joint_acc,'linewidth',2)
    yyaxis right
    plot(td_act(trial_num).opensim(:,10),'linewidth',2)
    hold on
    ylims = get(gca,'ylim');
    plot(repmat(td_act(trial_num).idx_goCueTime,2,1),ylims,'k--','linewidth',2)
    title 'Active shoulder flexion'
    set(gca,'box','off','tickdir','out')

    subplot 212
    yyaxis left
    %calculate acceleration
    joint_acc = gradient(td_pas(trial_num).opensim(:,3),1);
    plot(joint_acc,'linewidth',2)
    yyaxis right
    plot(td_pas(trial_num).opensim(:,10),'linewidth',2)
    hold on
    ylims = get(gca,'ylim');
    plot(repmat(td_pas(trial_num).idx_bumpTime,2,1),ylims,'k--','linewidth',2)
    title 'Passive shoulder flexion'
    set(gca,'box','off','tickdir','out')
end