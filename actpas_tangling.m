%% Preprocess trial data
    [~,td] = getTDidx(trial_data,'result','R');
    
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
    td = addFiringRates(td,struct('array','S1'));
    
    % trim active
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    % clean nans out...?
    nanners = isnan(cat(1,td_act.target_direction));
    td_act = td_act(~nanners);
    % td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',30});
    td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_endTime',0});
    td_act = binTD(td_act,5);

    % trim passive
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',30});
    td_pas = binTD(td_pas,5);

    % put active and passive back together
    td = cat(2,td_act,td_pas);

%% Create data struct for tangling code
    S1_tangle_struct = struct('A',{td_act.S1_FR});
    EMG_tangle_struct = struct('A',{td_act.emg});

    [Q_S1,out_S1] = tangleAnalysis(S1_tangle_struct,0.05,'softenNorm',1);
    [Q_EMG,out_EMG] = tangleAnalysis(EMG_tangle_struct,0.05,'softenNorm',0);

%% Scatter plot...
    figure
    scatter(Q_EMG,Q_S1)
    hold on
    set(gca,'box','off','tickdir','out')
    axis([0 inf 0 inf])
    % axis equal
