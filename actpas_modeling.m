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

%% Try fabricating trial_data with linear models based on muscles
    % get models for force and velocity from actpas data
    % opensim_idx = find(contains(td(1).opensim_names,'_moment'));
    % opensim_idx = find(contains(td(1).opensim_names,'_vel'));
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
    % opensim_idx = find(contains(td(1).opensim_names,'_vel') & ~contains(td(1).opensim_names,'wrist') & ~contains(td(1).opensim_names,'radial'));
    % opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
    % opensim_idx = find(contains(td(1).opensim_names,'_vel') | contains(td(1).opensim_names,'_moment'));
    % opensim_idx = find( (contains(td(1).opensim_names,'_vel') | contains(td(1).opensim_names,'_moment')) & ~contains(td(1).opensim_names,'wrist') & ~contains(td(1).opensim_names,'radial'));
    % get active and passive trial indices
    active_trials = getTDidx(td,'ctrHoldBump',false);
    passive_trials = getTDidx(td,'ctrHoldBump',true);

    [td,model_info_full] = getModel(td,struct('model_type','glm',...
        'model_name','S1_muscle','in_signals',...
        {{'opensim',opensim_idx}},...
        'out_signals',{'S1_FR'}));
    [td,model_info_act] = getModel(td,struct('model_type','glm',...
        'model_name','S1_muscle_act','in_signals',...
        {{'opensim',opensim_idx}},...
        'out_signals',{'S1_FR'},'train_idx',active_trials));
    [td,model_info_pas] = getModel(td,struct('model_type','glm',...
        'model_name','S1_muscle_pas','in_signals',...
        {{'opensim',opensim_idx}},...
        'out_signals',{'S1_FR'},'train_idx',passive_trials));

    [td_act,model_info_act] = getModel(td(active_trials),struct('model_type','glm',...
        'model_name','S1_muscle_combined','in_signals',...
        {{'opensim',opensim_idx}},...
        'out_signals',{'S1_FR'}));
    [td_pas,model_info_pas] = getModel(td(passive_trials),struct('model_type','glm',...
        'model_name','S1_muscle_combined','in_signals',...
        {{'opensim',opensim_idx}},...
        'out_signals',{'S1_FR'}));
    td = [td_act,td_pas];
