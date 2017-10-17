function [separability,sep_train] = test_sep(td,signals,num_dim)

% td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
[td,~] = getPCA(td,struct('signals',signals));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);

% clean nans out...?
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);

% plot active as filled, passive as open
bump_colors = linspecer(4);
if strfind(signals{1}{1},'_spikes')
    pca_name = [strtok(signals{1}{1},'_spikes') '_pca'];
else
    pca_name = [signals{1}{1} '_pca'];
end
S1_pca_act = cat(1,td_act.(pca_name));
S1_pca_pas = cat(1,td_pas.(pca_name));
act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;

% fix num_dim if few neurons
num_dim = min(num_dim,size(td(1).(pca_name),2));

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
mdl = fitcdiscr(S1_pca(train_idx,1:num_dim),actpas(train_idx));
class = predict(mdl,S1_pca(test_idx,1:num_dim));
separability = sum(class == actpas(test_idx))/sum(test_idx);

class_train = predict(mdl,S1_pca(train_idx,1:num_dim));
sep_train = sum(class_train == actpas(train_idx))/sum(train_idx);

w = mdl.Sigma\diff(mdl.Mu)';
figure
hold all
scatter(S1_pca_act(:,1:num_dim)*w,1:length(S1_pca_act),50,bump_colors(act_dir_idx,:),'filled')
scatter(S1_pca_pas(:,1:num_dim)*w,1:length(S1_pca_pas),50,bump_colors(pas_dir_idx,:))
set(gca,'box','off','tickdir','out')
