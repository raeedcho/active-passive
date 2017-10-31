%% extract OOR reaches for R script
[~,td] = getTDidx(trial_data,'result','R');

[fr,thv,thf] = obs_window_pds(td);

scatter(thf,thv,'o')
set(gca,'xlim',[-pi pi],'ylim',[-pi pi])

folder = 'C:\Users\rhc307\Projects\limblab\data-raeed\ForceKin\Han\20170203\OOR_CSV\';
filename = 'Han_20170203_OOR_reaches';
% csvwrite([folder filename '.csv'],[thv thf fr])

% plot tuning surfaces
plot_tuning_surfaces(fr,thv,thf);

% extract residuals into csv
[resid_lin, resid_non] = extract_resid(fr,thv,thf);
csvwrite([folder 'Han_20170203_OOR_resid_lin.csv'],[thv thf resid_lin])
csvwrite([folder 'Han_20170203_OOR_resid_nonlin.csv'],[thv thf resid_non])

%% try removing nonlinear force-velocity combination
% [~,td] = getTDidx(trial_data,'result','R');
% 
% [fr,thv,thf] = obs_window_pds(td);
% 
% % get linear models
% [lin_models] = get_lin_models(fr,thv,thf);
% 
% % fit active/passive stuff
% [~,td_actpas] = getTDidx(trial_data_actpas,'result','R');
% 
% % test separability of original actpas
% [sep_actual,sep_train_actual] = test_sep(td_actpas);
% 
% % predict td from models
% td_actpas_pred = predict_td_from_models(td_actpas,lin_models);
% 
% % test separability of lin model
% [sep_pred, sep_train_pred] = test_sep(td_actpas_pred);

%% test removal of residual with fine linear model
[~,td] = getTDidx(trial_data,'result','R');

% get linear models
[lin_models] = get_lin_models_from_td(td);

% fit active/passive stuff
[~,td_actpas] = getTDidx(trial_data_actpas,'result','R');

% test separability of original actpas
[sep_actual,sep_train_actual] = test_sep(td_actpas);

% predict td from models
td_actpas_pred = predict_td_from_models(td_actpas,lin_models);

% test separability of lin model
[sep_pred, sep_train_pred] = test_sep(td_actpas_pred);

%% check linear models
[~,td] = getTDidx(trial_data,'result','R');

[fr,thv,thf] = obs_window_pds(td);


% get linear models
[lin_models] = get_lin_models(fr,thv,thf);

plotDecomposedFR(fr,thv,thf,lin_models,true);

% predict OOR with lin models
% td_pred = predict_td_from_models(td,lin_models);
% 
% [fr,thv,thf] = obs_window_pds(td_pred);
% 
% plotDecomposedFR(fr,thv,thf,lin_models,false

%% Check linear model extraction
[~,td] = getTDidx(trial_data,'result','R');

% [fr,thv,thf] = obs_window_pds(td);

% get linear models
[lin_models] = get_lin_models_from_td(td);

% predict td from models
td_pred = predict_td_from_models(td,lin_models);
[fr,thv,thf] = obs_window_pds(td_pred);
[lin_models] = get_lin_models(fr,thv,thf);
plotDecomposedFR(fr,thv,thf,lin_models,true)

%% try linear model on just actpas
[~,td] = getTDidx(trial_data_actpas,'result','R');
td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime'));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_goCueTime',35},{'idx_goCueTime',50});
% td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
td_act = binTD(td_act,15);
% clean nans out...?
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);
tgt_dir = cat(1,td_act.target_direction);

[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_pas = binTD(td_pas,15);
bump_dir = cat(1,td_pas.bumpDir)*pi/180;

td = cat(2,td_act,td_pas);
thv = [tgt_dir;bump_dir];
thf = [-tgt_dir;bump_dir];
fr = cat(1,td.S1_spikes);
fr_pred = zeros(size(fr));

%fit and evaluate linear models
for neuron_ctr = 1:size(fr,2)
    dir_mat = [cos(thv) sin(thv) cos(thf) sin(thf)];
    fit_linear = LinearModel.fit(dir_mat,fr(:,neuron_ctr));

    fr_pred(:,neuron_ctr) = fit_linear.predict(dir_mat);
end
fr_pred(fr_pred<0) = 0;

% assign linear predictions
for trial_ctr = 1:length(td)
    td(trial_ctr).S1_spikes = fr_pred(trial_ctr,:);
end

% test
td = sqrtTransform(td,'S1_spikes');
[td,~] = getPCA(td,struct('signals',{'S1_spikes'}));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);

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
% S1_pca = cat(1,S1_pca_act,S1_pca_pas);
% actpas = [ones(length(S1_pca_act),1);zeros(length(S1_pca_pas),1)];
% [train_idx,test_idx] = crossvalind('LeaveMOut',length(actpas),floor(length(actpas)/10));
% mdl = fitcdiscr(S1_pca(train_idx,:),actpas(train_idx));
% class = predict(mdl,S1_pca(test_idx,:));
% separability = sum(class == actpas(test_idx))/sum(test_idx);
% 
% class_train = predict(mdl,S1_pca(train_idx,:));
% sep_train = sum(class_train == actpas(train_idx))/sum(train_idx);
% 
% w = mdl.Sigma\diff(mdl.Mu)';
% figure
% hold all
% scatter(S1_pca_act*w,1:length(S1_pca_act),50,bump_colors(act_dir_idx,:),'filled')
% scatter(S1_pca_pas*w,1:length(S1_pca_pas),50,bump_colors(pas_dir_idx,:))
% set(gca,'box','off','tickdir','out')
