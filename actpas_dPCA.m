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

%% try dPCA
    td = smoothSignals(td,struct('signals',{{'S1_spikes'}},'calc_rate',true,'kernel_SD',0.05));
    [td,dPCA_info] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'S1_spikes'}},'num_dims',8,'do_plot',true,...
                                                                        'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                        'dpca_plot_fcn',@dpca_plot_actpas));
    
%% Visualize direction and active/passive dimensions
    active_trials = getTDidx(td,'ctrHoldBump',false);
    passive_trials = getTDidx(td,'ctrHoldBump',true);
    td_act_avg = trialAverage(td(active_trials),'target_direction');
    td_pas_avg = trialAverage(td(passive_trials),'bumpDir');

    % td_act_avg = td(active_trials);
    % td_pas_avg = td(passive_trials);

    % signal = 'S1_spikes';
    % signal = 'glm_S1_muscle';
    signal = 'opensim_pca';
    % W = dPCA_info.W(:,[2 5 3]);
    % W = musc_dpca.W(:,[1 4 3]);
    W = muscVel_dPCA.W(:,[1 4 3])*300;
    mean_spikes = mean(cat(1,td.(signal)));
    dir_colors = linspecer(4);
    num_bins = 68;
    fignum = figure;
    % animation stuff
    scatter(randn(1,4),randn(1,4),[],dir_colors,'filled')
    f = getframe(fignum);
    [im,map] = rgb2ind(f.cdata,256,'nodither');
    im(1,1,1,num_bins) = 0;
    % set(gca,'visible','off')
    for time = 1:num_bins
        clf
        for trial = 1:length(td_pas_avg)
            X = td_pas_avg(trial).(signal);
            Xcen = X-mean_spikes;
            Z = Xcen * W;
    
            dir_idx = td_pas_avg(trial).bumpDir/90+1;

            subplot(1,2,1)
            plot(Z(1:time,1),Z(1:time,2),'--','linewidth',2,'color',dir_colors(dir_idx,:))
            hold on
            plot(Z(time,1),Z(time,2),'x','markersize',10,'color',dir_colors(dir_idx,:))

            subplot(1,2,2)
            plot(Z(1:time,1),Z(1:time,3),'--','linewidth',2,'color',dir_colors(dir_idx,:))
            hold on
            plot(Z(time,1),Z(time,3),'x','markersize',10,'color',dir_colors(dir_idx,:))
        end
        for trial = 1:length(td_act_avg)
            X = td_act_avg(trial).(signal);
            Xcen = X-mean_spikes;
            Z = Xcen * W;
    
            dir_idx = td_act_avg(trial).target_direction/(pi/2)+1;
            dir_idx = round(dir_idx);

            subplot(1,2,1)
            plot(Z(1:time,1),Z(1:time,2),'-','linewidth',2,'color',dir_colors(dir_idx,:))
            hold on
            plot(Z(time,1),Z(time,2),'.','markersize',40,'color',dir_colors(dir_idx,:))

            subplot(1,2,2)
            plot(Z(1:time,1),Z(1:time,3),'-','linewidth',2,'color',dir_colors(dir_idx,:))
            hold on
            plot(Z(time,1),Z(time,3),'.','markersize',40,'color',dir_colors(dir_idx,:))
        end
        subplot(1,2,1)
        axis equal
        axis([-30 30 -30 30 -30 30])
        set(gca,'box','off','tickdir','out')
        subplot(1,2,2)
        axis equal
        axis([-30 30 -30 30 -30 30])
        set(gca,'box','off','tickdir','out')

        % get animation pieces
        f = getframe(fignum);
        im(:,:,1,time) = rgb2ind(f.cdata,map,'nodither');

        pause(0.1)
    end
    imwrite(im,map,'test.gif','DelayTime',0,'LoopCount',0)

%% Try fabricating trial_data with linear models based on handle kinematics and force
    % get models for force and velocity from actpas data
    [td,model_info] = getModel(td,struct('model_type','linmodel',...
        'model_name','S1_handle','in_signals',{{'vel';'force'}},...
        'out_signals',{'S1_spikes'}));
    
    % [td,dPCA_info] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'linmodel_S1_handle'}},'num_dims',4,'do_plot',true,...
    %                                                                     'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
    %                                                                     'dpca_plot_fcn',@dpca_plot_actpas));

%% Try fabricating trial_data with linear models based on muscles
    % get covariates to base simulated neurons on
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
    td = getPCA(td,struct('signals',{{'opensim',opensim_idx}},'do_plot',false));
    % get dpca for muscle velocity signals
    [~,muscVel_dPCA] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'opensim_pca'}},'num_dims',8,'do_plot',true,...
                                                                        'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                        'dpca_plot_fcn',@dpca_plot_actpas));

    % get active and passive trial indices
    active_trials = getTDidx(td,'ctrHoldBump',false);
    passive_trials = getTDidx(td,'ctrHoldBump',true);

    % % linear models
    % [td,model_info_full] = getModel(td,struct('model_type','linmodel',...
    %     'model_name','S1_muscle','in_signals',...
    %     {{'opensim_pca',1:6}},...
    %     'out_signals',{'S1_FR'}));
    % [td,model_info_act] = getModel(td,struct('model_type','linmodel',...
    %     'model_name','S1_muscle_act','in_signals',...
    %     {{'opensim_pca',1:6}},...
    %     'out_signals',{'S1_FR'},'train_idx',active_trials));
    % [td,model_info_pas] = getModel(td,struct('model_type','linmodel',...
    %     'model_name','S1_muscle_pas','in_signals',...
    %     {{'opensim_pca',1:6}},...
    %     'out_signals',{'S1_FR'},'train_idx',passive_trials));
    % 
    % % get dpca for muscle-based neurons
    % getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'linmodel_S1_muscle'}},'num_dims',8,'do_plot',true,...
    %                                                                     'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
    %                                                                     'dpca_plot_fcn',@dpca_plot_actpas));
    % % project into original dPCA space
    % getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'linmodel_S1_muscle'}},'num_dims',8,'do_plot',true,...
    %                                                                     'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
    %                                                                     'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));
    % % get dpca for active-trained muscle-based neurons, projected into original space
    % getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'linmodel_S1_muscle_act'}},'num_dims',8,'do_plot',true,...
    %                                                                     'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
    %                                                                     'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));
    % % get dpca for passive-trained muscle-based neurons, projected into original space
    % getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'linmodel_S1_muscle_pas'}},'num_dims',8,'do_plot',true,...
    %                                                                     'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
    %                                                                     'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));

    % glms
    [td,model_info_full] = getModel(td,struct('model_type','glm',...
        'model_name','S1_muscle','in_signals',...
        {{'opensim_pca',1:6}},...
        'out_signals',{'S1_FR'}));
    [td,model_info_act] = getModel(td,struct('model_type','glm',...
        'model_name','S1_muscle_act','in_signals',...
        {{'opensim_pca',1:6}},...
        'out_signals',{'S1_FR'},'train_idx',active_trials));
    [td,model_info_pas] = getModel(td,struct('model_type','glm',...
        'model_name','S1_muscle_pas','in_signals',...
        {{'opensim_pca',1:6}},...
        'out_signals',{'S1_FR'},'train_idx',passive_trials));
    [td_act,model_info_act] = getModel(td(active_trials),struct('model_type','glm',...
        'model_name','S1_muscle_combined','in_signals',...
        {{'opensim_pca',1:6}},...
        'out_signals',{'S1_FR'}));
    [td_pas,model_info_pas] = getModel(td(passive_trials),struct('model_type','glm',...
        'model_name','S1_muscle_combined','in_signals',...
        {{'opensim_pca',1:6}},...
        'out_signals',{'S1_FR'}));
    td = [td_act,td_pas];
    
    % get dpca for muscle-based neurons
    [~,musc_dpca] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'glm_S1_muscle'}},'num_dims',8,'do_plot',true,...
                                                                        'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                        'dpca_plot_fcn',@dpca_plot_actpas));
    % project into original dPCA space
    getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'glm_S1_muscle'}},'num_dims',8,'do_plot',true,...
                                                                        'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                        'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));
    % get dpca for active-trained muscle-based neurons, projected into original space
    getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'glm_S1_muscle_act'}},'num_dims',8,'do_plot',true,...
                                                                        'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                        'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));
    % get dpca for passive-trained muscle-based neurons, projected into original space
    getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'glm_S1_muscle_pas'}},'num_dims',8,'do_plot',true,...
                                                                        'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                        'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));
    % get dpca for combined active-trained and passive-trained muscle-based neurons, projected into original space
    getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'glm_S1_muscle_combined'}},'num_dims',8,'do_plot',true,...
                                                                        'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                        'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));

%% Plot average neural firing rates against muscle-based predictions
    % split td and average
    td_avg = trialAverage(td,{'ctrHoldBump','target_direction'});
    [~,td_act] = getTDidx(td_avg,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td_pas,'ctrHoldBump',true);

    % plot for each neuron
    dir_colors = linspecer(4);
    for neuron_idx = 1:size(td(1).S1_spikes,2)
        h=figure;
        for dir_idx = 1:length(td_act)
            subplot(4,2,(dir_idx-1)*2+1)
            plot(td_act(dir_idx).S1_spikes(:,neuron_idx),'-','linewidth',2,'color',dir_colors(dir_idx,:))
            hold on
            plot(td_act(dir_idx).glm_S1_muscle(:,neuron_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
            axis tight
            set(gca,'box','off','tickdir','out','ylim',[0 inf])
            ylim = get(gca,'ylim');

            subplot(4,2,(dir_idx-1)*2+2)
            plot(td_pas(dir_idx).S1_spikes(:,neuron_idx),'-','linewidth',2,'color',dir_colors(dir_idx,:))
            hold on
            plot(td_pas(dir_idx).glm_S1_muscle(:,neuron_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
            axis tight
            set(gca,'box','off','tickdir','out','ylim',ylim)
        end
        waitfor(h)
    end
