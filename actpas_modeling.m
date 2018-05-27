%% setup
    [~,td] = getTDidx(trial_data,'result','R');
    
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
    
    % get still handle data (no word for start of center hold)
    minCH = min(cat(1,td.ctrHold));
    bin_size = td(1).bin_size;
    still_bins = floor(minCH/bin_size);
    
    % Get td_act and td_pas
    num_bins_before = floor(still_bins/2);
    num_bins_after = 30;
    
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

    % Add new firing rates signal
    td = addFiringRates(td,struct('array','S1'));
    % smooth signals at 50 ms
    % td = sqrtTransform(td,'S1_spikes');
    % td = smoothSignals(td,struct('signals',{{'S1_spikes'}},'calc_rate',true,'kernel_SD',0.05));

    % clean up
    clearvars -except trial_data td

%% Examine efficacy of a hybrid active/passive model against individual models
    % get models for force and velocity from actpas data
    % opensim_idx = find(contains(td(1).opensim_names,'_moment'));
    % opensim_idx = find(contains(td(1).opensim_names,'_vel'));
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
    td = getPCA(td,struct('signals',{{'opensim',opensim_idx}},'do_plot',true));
    % opensim_idx = find(contains(td(1).opensim_names,'_vel') & ~contains(td(1).opensim_names,'wrist') & ~contains(td(1).opensim_names,'radial'));
    % opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
    % opensim_idx = find(contains(td(1).opensim_names,'_vel') | contains(td(1).opensim_names,'_moment'));
    % opensim_idx = find( (contains(td(1).opensim_names,'_vel') | contains(td(1).opensim_names,'_moment')) & ~contains(td(1).opensim_names,'wrist') & ~contains(td(1).opensim_names,'radial'));
    % get active and passive trial indices
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    [train_act,test_act] = crossvalind('HoldOut',length(td_act),0.5);
    [train_pas,test_pas] = crossvalind('HoldOut',length(td_pas),0.5);
    td_train_act = td_act(train_act);
    td_train_pas = td_pas(train_pas);
    td_test_act = td_act(test_act);
    td_test_pas = td_pas(test_pas);
    td_train = [td_train_act td_train_pas];
    td_test = [td_test_act td_test_pas];

    % first try training on active and testing on active and passive
    [td_train_act,act_model] = getModel(td_train_act,struct('model_type','glm',...
        'model_name','S1_act','in_signals',...
        {{'opensim_pca',1:6}},...
        'out_signals',{'S1_FR'}));
    act_model.eval_metric = 'pr2';
    act_model.num_boots = 1;
    td_test_act = getModel(td_test_act,act_model);
    act_within_pR2 = evalModel(td_test_act,act_model);
    td_test_pas = getModel(td_test_pas,act_model);
    act_across_pR2 = evalModel(td_test_pas,act_model);

    % then train on passive
    [td_train_pas,pas_model] = getModel(td_train_pas,struct('model_type','glm',...
        'model_name','S1_pas','in_signals',...
        {{'opensim_pca',1:6}},...
        'out_signals',{'S1_FR'}));
    pas_model.eval_metric = 'pr2';
    pas_model.num_boots = 1;
    td_test_act = getModel(td_test_act,pas_model);
    pas_across_pR2 = evalModel(td_test_act,pas_model);
    td_test_pas = getModel(td_test_pas,pas_model);
    pas_within_pR2 = evalModel(td_test_pas,pas_model);

    % then some sort of hybrid...
    [td_train,hybrid_model] = getModel(td_train,struct('model_type','glm',...
        'model_name','S1_hybrid','in_signals',...
        {{'opensim_pca',1:6}},...
        'out_signals',{'S1_FR'}));
    hybrid_model.eval_metric = 'pr2';
    hybrid_model.num_boots = 1;
    td_test_act = getModel(td_test_act,hybrid_model);
    hybrid_act_pR2 = evalModel(td_test_act,hybrid_model);
    td_test_pas = getModel(td_test_pas,hybrid_model);
    hybrid_pas_pR2 = evalModel(td_test_pas,hybrid_model);

    % plot them against each other
    figure
    scatter(act_within_pR2,pas_within_pR2,'filled')
    hold on
    plot([-0.05 0.2],[-0.05 0.2],'k--','linewidth',2)
    plot([0 0],[-0.05 0.2],'k-','linewidth',2)
    plot([-0.05 0.2],[0 0],'k-','linewidth',2)
    set(gca,'box','off','tickdir','out','xlim',[-0.05 0.2],'ylim',[-0.05 0.2])
    axis equal
    xlabel 'Active pR^2'
    ylabel 'Passive pR^2'
    % figure
    % scatter(act_within_pR2,act_across_pR2,'filled')
    % hold on
    % plot([-0.05 0.2],[-0.05 0.2],'k--','linewidth',2)
    % plot([0 0],[-0.05 0.2],'k-','linewidth',2)
    % plot([-0.05 0.2],[0 0],'k-','linewidth',2)
    % set(gca,'box','off','tickdir','out','xlim',[-0.05 0.2],'ylim',[-0.05 0.2])
    % figure
    % scatter(pas_within_pR2,pas_across_pR2,'filled')
    % hold on
    % plot([-0.05 0.2],[-0.05 0.2],'k--','linewidth',2)
    % plot([0 0],[-0.05 0.2],'k-','linewidth',2)
    % plot([-0.05 0.2],[0 0],'k-','linewidth',2)
    % set(gca,'box','off','tickdir','out','xlim',[-0.05 0.2],'ylim',[-0.05 0.2])
    figure
    scatter(hybrid_act_pR2,hybrid_pas_pR2,'filled')
    hold on
    plot([-0.05 0.2],[-0.05 0.2],'k--','linewidth',2)
    plot([0 0],[-0.05 0.2],'k-','linewidth',2)
    plot([-0.05 0.2],[0 0],'k-','linewidth',2)
    set(gca,'box','off','tickdir','out','xlim',[-0.05 0.2],'ylim',[-0.05 0.2])
    axis equal
    xlabel 'Active pR^2'
    ylabel 'Passive pR^2'
