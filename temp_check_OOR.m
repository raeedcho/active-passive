%% extract OOR reaches for R script
[~,td] = getTDidx(trial_data_OOR,'result','R');

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
[~,td] = getTDidx(trial_data_OOR,'result','R');

[fr,thv,thf] = obs_window_pds(td);

% get linear models
[lin_models] = get_lin_models(fr,thv,thf);

% fit active/passive stuff
[~,td_actpas] = getTDidx(trial_data_actpas,'result','R');

% test separability of original actpas
[sep_actual,sep_train_actual] = test_sep(td_actpas);

% predict td from models
td_actpas_pred = predict_td_from_models(td_actpas,lin_models);

% test separability of lin model
[sep_pred, sep_train_pred] = test_sep(td_actpas_pred);

%% check linear models
[~,td] = getTDidx(trial_data_OOR,'result','R');

[fr,thv,thf] = obs_window_pds(td);


% get linear models
[lin_models] = get_lin_models(fr,thv,thf);

plotDecomposedFR(fr,thv,thf,lin_models,true);

% predict OOR with lin models
% td_pred = predict_td_from_models(td,lin_models);
% 
% [fr,thv,thf] = obs_window_pds(td_pred);
% 
% plotDecomposedFR(fr,thv,thf,lin_models,false);