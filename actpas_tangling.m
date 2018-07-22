%% Preprocess trial data
    [~,td] = getTDidx(trial_data,'result','R');
    
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

    % smooth signals
    td = smoothSignals(td,struct('signals','S1_spikes','kernel_SD',0.02,'calc_rate',true));
    td = smoothSignals(td,struct('signals','emg','kernel_SD',0.02,'calc_rate',false));
    
    % trim active
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    % clean nans out...?
    nanners = isnan(cat(1,td_act.target_direction));
    td_act = td_act(~nanners);
    % td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',30});
    td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_endTime',0});
    % td_act = binTD(td_act,5);

    % trim passive
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',30});
    % td_pas = binTD(td_pas,5);

    % put active and passive back together
    td = cat(2,td_act,td_pas);

%% Sanitize data based on kinematics
    % check time to completion for active movements
    duration_list = zeros(length(td_act),1);
    figure
    for i = 1:length(td_act)
        duration_list(i) = size(td_act(i).pos,1);
        if duration_list(i)<45 % arbitrarily chosen threshold based on this plot
            plot(td_act(i).pos(:,1),td_act(i).pos(:,2),'-k')
        else
            plot(td_act(i).pos(:,1),td_act(i).pos(:,2),'-r')
        end
        hold on
    end
    axis equal
    set(gca,'box','off','tickdir','out')

    % remove long duration traces
    td_act(duration_list>=45) = [];

%% Trial average signals
    td_act = trialAverage(td_act,'target_direction',struct('do_stretch',true,'num_samp',31));
    td_pas = trialAverage(td_pas,'bumpDir');

%% Create data struct for tangling code
    S1_tangle_struct = struct('A',{td_act.S1_spikes});
    EMG_tangle_struct = struct('A',{td_act.emg});
    [Q_S1_act,out_S1_act] = tangleAnalysis(S1_tangle_struct,0.01,'softenNorm',5,'timeStep',1);
    [Q_EMG_act,out_EMG_act] = tangleAnalysis(EMG_tangle_struct,0.01,'softenNorm',0,'timeStep',1);

    S1_tangle_struct = struct('A',{td_pas.S1_spikes});
    EMG_tangle_struct = struct('A',{td_pas.emg});
    [Q_S1_pas,out_S1_pas] = tangleAnalysis(S1_tangle_struct,0.01,'softenNorm',5,'timeStep',1);
    [Q_EMG_pas,out_EMG_pas] = tangleAnalysis(EMG_tangle_struct,0.01,'softenNorm',0,'timeStep',1);


%% Scatter plot...
    lims = max([Q_S1_pas;Q_S1_act]);
    figure
    scatter(Q_S1_pas,Q_S1_act,[],'k','filled')
    hold on
    plot([0 lims],[0 lims],'--k','linewidth',2)
    set(gca,'box','off','tickdir','out')
    axis equal
    axis([0 lims 0 lims])
    xlabel 'Tangling during passive'
    ylabel 'Tangling during active'

%% motor cortex...
    td_motor = load('~/Projects/limblab/data-td/MattData/Chewie_CO_FF_2016-10-07.mat','trial_data');
    td_motor = td_motor.trial_data;
    [~,td_motor_bl] = getTDidx(td_motor,'result','R','epoch','BL');
    td_motor_bl = getMoveOnsetAndPeak(td_motor_bl,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

    % smooth signals
    td_motor_bl = smoothSignals(td_motor_bl,struct('signals','M1_spikes','kernel_SD',0.02,'calc_rate',true));
    td_motor_bl = smoothSignals(td_motor_bl,struct('signals','PMd_spikes','kernel_SD',0.02,'calc_rate',true));
    % trim signals
    td_motor_bl = trimTD(td_motor_bl,{'idx_movement_on',-30},{'idx_movement_on',30});

    % check time to completion for active movements
    dist_list = zeros(length(td_motor_bl),1);
    figure
    for i = 1:length(td_motor_bl)
        dist_list(i) = sqrt(sum((td_motor_bl(i).pos(end,:)-td_motor_bl(i).pos(1,:)).^2,2));
        if dist_list(i)<3 % arbitrarily chosen threshold based on this plot
            plot(td_motor_bl(i).pos(:,1),td_motor_bl(i).pos(:,2),'-r')
        else
            plot(td_motor_bl(i).pos(:,1),td_motor_bl(i).pos(:,2),'-k')
        end
        hold on
    end
    axis equal
    set(gca,'box','off','tickdir','out')

    td_motor_bl(dist_list<3) = [];
    
    % trial average and tangling
    td_motor_bl = trialAverage(td_motor_bl,'target_direction');
    M1_tangle_struct = struct('A',{td_motor_bl.M1_spikes});
    PMd_tangle_struct = struct('A',{td_motor_bl.PMd_spikes});

    [Q_M1,out_M1] = tangleAnalysis(M1_tangle_struct,0.01,'softenNorm',5,'timeStep',1);
    [Q_PMd,out_PMd] = tangleAnalysis(PMd_tangle_struct,0.01,'softenNorm',5,'timeStep',1);

    lims = max([Q_M1;Q_PMd]);
    figure
    scatter(Q_M1,Q_PMd,[],'k','filled')
    hold on
    plot([0 lims],[0 lims],'--k','linewidth',2)
    set(gca,'box','off','tickdir','out')
    axis equal
    axis([0 lims 0 lims])
    xlabel 'M1 Tangle'
    ylabel 'PMd Tangle'
