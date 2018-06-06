%% setup
    [~,td] = getTDidx(trial_data,'result','R');
    
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
    
    % Get td_act and td_pas
    num_bins_before = 0;
    num_bins_after = 15;
    
    % prep td_act
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    td_act = trimTD(td_act,{'idx_movement_on',-num_bins_before},{'idx_movement_on',num_bins_after});
    % clean nans out...?
    nanners = isnan(cat(1,td_act.target_direction));
    td_act = td_act(~nanners);
    
    % prep td_pas
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',-num_bins_before},{'idx_bumpTime',num_bins_after});
    % move bumpDir into target_direction for passive td
    if floor(td_pas(1).bumpDir) == td_pas(1).bumpDir
        % probably in degrees
        multiplier = pi/180;
    else
        warning('bumpDir may be in radians')
        multiplier = 1;
    end
    for trial = 1:length(td_pas)
        td_pas(trial).target_direction = td_pas(trial).bumpDir*multiplier;
    end
    
    % concatenate tds
    td = cat(2,td_act,td_pas);

    % smooth signals at 50 ms and calculate firing rate by dividing by bin size
    td = sqrtTransform(td,'S1_spikes');
    td = smoothSignals(td,struct('signals',{{'S1_spikes'}},'calc_rate',true,'kernel_SD',0.05));

    % clean up
    clearvars -except trial_data td

    % set up colormap
    cm_viridis = viridis(200);

%% get pairwise corr in the two epochs
    % get active and passive trial indices
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

    [rho_act,sort_idx] = pairwiseCorr(td_act,struct('signals',{{'S1_spikes'}},'cluster_order',true));
    [rho_pas] = pairwiseCorr(td_pas,struct('signals',{{'S1_spikes'}}));

    figure
    subplot(2,1,1)
    imagesc(rho_act(sort_idx,sort_idx))
    clim = get(gca,'clim');
    colorbar
    colormap(cm_viridis);
    axis square
    title 'S1 neural correlation - Active'
    subplot(2,1,2)
    imagesc(rho_pas(sort_idx,sort_idx),clim)
    colorbar
    colormap(cm_viridis);
    axis square
    title 'S1 neural correlation - Passive'

%% correlation analysis in muscle space
    % get covariates to base simulated neurons on
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));

    % get active and passive trial indices
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

    [rho_act,sort_idx] = pairwiseCorr(td_act,struct('signals',{{'opensim',opensim_idx}},'cluster_order',true));
    [rho_pas] = pairwiseCorr(td_pas,struct('signals',{{'opensim',opensim_idx}}));

    figure
    subplot(2,1,1)
    imagesc(rho_act(sort_idx,sort_idx))
    clim = get(gca,'clim');
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Muscle velocity correlation - Active'
    subplot(2,1,2)
    imagesc(rho_pas(sort_idx,sort_idx),clim)
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Muscle velocity correlation - Passive'

%% correlation analysis in joint space
    % get covariates to base simulated neurons on
    opensim_idx = find(contains(td(1).opensim_names,'_vel'));

    % get active and passive trial indices
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

    [rho_act,sort_idx] = pairwiseCorr(td_act,struct('signals',{{'opensim',opensim_idx}},'cluster_order',true));
    [rho_pas] = pairwiseCorr(td_pas,struct('signals',{{'opensim',opensim_idx}}));

    figure
    subplot(2,1,1)
    imagesc(rho_act(sort_idx,sort_idx))
    clim = get(gca,'clim');
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Joint velocity correlation - Active'
    subplot(2,1,2)
    imagesc(rho_pas(sort_idx,sort_idx),clim)
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Joint velocity correlation - Passive'

%% correlation analysis of muscle based neurons
    % get covariates to base simulated neurons on
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
    td = getPCA(td,struct('signals',{{'opensim',opensim_idx}},'do_plot',false));

    % fit model
    [td,model_info_full] = getModel(td,struct('model_type','glm',...
        'model_name','S1_muscle','in_signals',...
        {{'opensim_pca',1:6}},...
        'out_signals',{'S1_spikes'}));

    % get active and passive trial indices
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

    [rho_act,sort_idx] = pairwiseCorr(td_act,struct('signals',{{'glm_S1_muscle'}},'cluster_order',true));
    [rho_pas] = pairwiseCorr(td_pas,struct('signals',{{'glm_S1_muscle'}}));

    figure
    subplot(2,1,1)
    imagesc(rho_act(sort_idx,sort_idx))
    clim = get(gca,'clim');
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Muscle-based neural correlation - Active'
    subplot(2,1,2)
    imagesc(rho_pas(sort_idx,sort_idx),clim)
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Muscle-based neural correlation - Passive'
