%% extract OOR reaches for R script
[~,td] = getTDidx(trial_data,'result','R');

[fr,thv,thf] = obs_window_pds(td,struct('signals',{{'S1_spikes'}}));

% tesselate...
% get linear models
[lin_models] = get_lin_models(fr,thv,thf);

plotDecomposedFR(fr,thv,thf,lin_models,true);


