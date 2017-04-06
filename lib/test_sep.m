function [separability,sep_train] = test_sep(td)

% get only reward trials
[~,td] = getTDidx(td,'result','R');

% get movement onset
td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime'));

% extract relevant times
[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_goCueTime',35},{'idx_goCueTime',50});
% td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
td_act = binTD(td_act,15);

[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_pas = binTD(td_pas,15);
td = cat(2,td_act,td_pas);

% td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
td = sqrtTransform(td,'S1_spikes');
[td,~] = getPCA(td,struct('signals',{'S1_spikes'}));

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

% Find total separability
S1_pca = cat(1,S1_pca_act,S1_pca_pas);
actpas = [ones(length(S1_pca_act),1);zeros(length(S1_pca_pas),1)];
[train_idx,test_idx] = crossvalind('LeaveMOut',length(actpas),floor(length(actpas)/10));
mdl = fitcdiscr(S1_pca(train_idx,:),actpas(train_idx));
class = predict(mdl,S1_pca(test_idx,:));
separability = sum(class == actpas(test_idx))/sum(test_idx);

class_train = predict(mdl,S1_pca(train_idx,:));
sep_train = sum(class_train == actpas(train_idx))/sum(train_idx);

w = mdl.Sigma\diff(mdl.Mu)';
figure
hold all
scatter(S1_pca_act*w,1:length(S1_pca_act),50,bump_colors(act_dir_idx,:),'filled')
scatter(S1_pca_pas*w,1:length(S1_pca_pas),50,bump_colors(pas_dir_idx,:))