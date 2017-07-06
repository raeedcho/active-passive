
%% Look at forces and speeds throughout movement
[~,td] = getTDidx(trial_data_actpas,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false,'target_direction',pi/2);
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',12});

[~,td_pas] = getTDidx(td,'ctrHoldBump',true,'bumpDir',90);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',12});

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

figure
for trial_num = 1:length(td_act)
    
    %calculate acceleration
    act_speed = sqrt(sum(td_act(trial_num).vel.^2,2));
    act_force = sqrt(sum(td_act(trial_num).force.^2,2));
    
    subplot 211
    
    yyaxis left
    plot(act_speed,'-','linewidth',1)
    hold on
    
    yyaxis right
    plot(act_force,'-','linewidth',1)
    hold on
    
    ylims = get(gca,'ylim');
%     plot(repmat(td_act(trial_num).idx_goCueTime,2,1),ylims,'k--','linewidth',1)
    plot(repmat(td_act(trial_num).idx_movement_on,2,1),ylims,'r--','linewidth',1)
    plot(repmat(td_act(trial_num).idx_peak_speed,2,1),ylims,'g--','linewidth',1)
    
    
    title 'Active movements'
    set(gca,'box','off','tickdir','out')
end

for trial_num = 1:length(td_pas)
    %calculate acceleration
    pas_speed = sqrt(sum(td_pas(trial_num).vel.^2,2));
    pas_force = sqrt(sum(td_pas(trial_num).force.^2,2));
    
    subplot 212
    
    yyaxis left
    plot(pas_speed,'-','linewidth',1)
    hold on
    
    yyaxis right
    plot(pas_force,'-','linewidth',1)
    hold on
    
    ylims = get(gca,'ylim');
    plot(repmat(td_pas(trial_num).idx_bumpTime,2,1),ylims,'k--','linewidth',1)
    
    title 'Passive movements'
    set(gca,'box','off','tickdir','out')
end



%% Get PCA for act vs pas

test_sep(trial_data_actpas)

%% Try fabricating trial_data with linear models based on behavior
[~,td] = getTDidx(trial_data_actpas,'result','R');

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

% td = sqrtTransform(td,'S1_spikes');
[td,~] = getPCA(td,struct('signals',{'linmodel_S1'}));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);

% clean nans out...?
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);

% plot active as filled, passive as open
bump_colors = linspecer(4);
S1_pca_act = cat(1,td_act.linmodel_S1_pca);
S1_pca_pas = cat(1,td_pas.linmodel_S1_pca);
act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;

% plot dots
figure
hold all
scatter3(S1_pca_act(:,1),S1_pca_act(:,2),S1_pca_act(:,3),50,bump_colors(act_dir_idx,:),'filled')
scatter3(S1_pca_pas(:,1),S1_pca_pas(:,2),S1_pca_pas(:,3),50,bump_colors(pas_dir_idx,:))
axis equal

% Find total separability
S1_pca = cat(1,S1_pca_act,S1_pca_pas);
actpas = [ones(length(S1_pca_act),1);zeros(length(S1_pca_pas),1)];
[train_idx,test_idx] = crossvalind('LeaveMOut',length(actpas),floor(length(actpas)/10));
mdl = fitcdiscr(S1_pca(train_idx,1:4),actpas(train_idx)); % with linear model, only 4 linearly independent columns
class = predict(mdl,S1_pca(test_idx,1:4));
separability = sum(class == actpas(test_idx))/sum(test_idx);

class_train = predict(mdl,S1_pca(train_idx,1:4));
sep_train = sum(class_train == actpas(train_idx))/sum(train_idx);

w = mdl.Sigma\diff(mdl.Mu)';
figure
hold all
scatter(S1_pca_act(:,1:4)*w,1:length(S1_pca_act),50,bump_colors(act_dir_idx,:),'filled')
scatter(S1_pca_pas(:,1:4)*w,1:length(S1_pca_pas),50,bump_colors(pas_dir_idx,:))
set(gca,'box','off','tickdir','out')

%% get PCA traces
[~,td] = getTDidx(trial_data_actpas,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});

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

