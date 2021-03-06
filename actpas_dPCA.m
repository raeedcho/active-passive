%% setup
    [~,td] = getTDidx(trial_data,'result','R');
    
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
    
    % get still handle data (no word for start of center hold)
    minCH = min(cat(1,td.ctrHold));
    bin_size = td(1).bin_size;
    still_bins = floor(minCH/bin_size);
    
    % Get td_act and td_pas
    % num_bins_before = floor(still_bins/2);
    num_bins_before = 30;
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

    % bin tds
    % td = binTD(td,5);

    % Add new firing rates signal
    td = addFiringRates(td,struct('array','S1'));
    % smooth signals at 50 ms
    % td = sqrtTransform(td,'S1_spikes');
    td = smoothSignals(td,struct('signals',{{'S1_spikes'}},'calc_rate',true,'kernel_SD',0.05));

    % clean up
    clearvars -except trial_data td

    % colormap
    cm_viridis = viridis(200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at dPCA
    %% try dPCA on actual neurons
    [td,dPCA_info] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'S1_spikes'}},'num_dims',6,'do_plot',true,...
                                                                        'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                        'dpca_plot_fcn',@dpca_plot_actpas));
    
    %% Modeling neurons based on muscle or joint
        % % This section models neurons on behavior to try explaining active and passive activity
        % % test out modeling methods by going from joint to muscle
        % % spoiler...it works really well
        % % get covariates to base simulated neurons on
        % musc_idx = find(contains(td(1).opensim_names,'_muscVel'));
        % joint_idx = find(contains(td(1).opensim_names,'_vel'));
        % td = getPCA(td,struct('signals',{{'opensim',musc_idx}},'do_plot',false));
        %
        % % model type for later use...
        % model_type = 'nn';
        %
        % % models
        % td = getModel(td,struct('model_type',model_type,...
        %     'model_name','joint_musc','in_signals',...
        %     {{'opensim',joint_idx}},...
        %     'out_signals',{{'opensim_pca',1:6}}));
        %
        % % split td and average
        % td_avg = trialAverage(td,{'ctrHoldBump','target_direction'});
        % [~,td_act] = getTDidx(td_avg,'ctrHoldBump',false);
        % [~,td_pas] = getTDidx(td_avg,'ctrHoldBump',true);

        % % plot for each neuron
        % dir_colors = linspecer(4);
        % for musc_idx = 1:6
        %     h=figure;
        %     ax = zeros(2,1);
        %     for dir_idx = 1:length(td_act)
        %         ax(1) = subplot(4,2,(dir_idx-1)*2+1);
        %         plot(td_act(dir_idx).opensim_pca(:,musc_idx),'-','linewidth',2,'color',dir_colors(dir_idx,:))
        %         hold on
        %         plot(td_act(dir_idx).nn_joint_musc(:,musc_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
        %         % plot(td_act(dir_idx).S1_FR(:,neuron_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
        %         axis tight
        %         set(gca,'box','off','tickdir','out')

        %         ax(2) = subplot(4,2,(dir_idx-1)*2+2);
        %         plot(td_pas(dir_idx).opensim_pca(:,musc_idx),'-','linewidth',2,'color',dir_colors(dir_idx,:))
        %         hold on
        %         plot(td_pas(dir_idx).nn_joint_musc(:,musc_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
        %         % plot(td_pas(dir_idx).S1_FR(:,neuron_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
        %         axis tight
        %         set(gca,'box','off','tickdir','out')
        %     end
        %     linkaxes(ax,'xy')
        %     waitfor(h)
        % end

    %% Try fabricating trial_data with linear models based on muscles
        % get covariates to base simulated neurons on
        joint_idx = find(contains(td(1).opensim_names,'_vel'));
        % opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
        % td = getPCA(td,struct('signals',{{'opensim',opensim_idx}},'do_plot',false));
        % get dpca for muscle velocity signals
        [~,muscVel_dPCA] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{'opensim_pca'}},'num_dims',8,'do_plot',true,...
                                                                            'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                            'dpca_plot_fcn',@dpca_plot_actpas));

        % get active and passive trial indices
        active_trials = getTDidx(td,'ctrHoldBump',false);
        passive_trials = getTDidx(td,'ctrHoldBump',true);

        % model type for later use...
        model_type = 'nn';

        % glms
        % [td,model_info_full] = getModel(td,struct('model_type',model_type,...
        %     'model_name','S1_muscle','in_signals',...
        %     {{'opensim',joint_idx}},...
        %     'out_signals',{'S1_FR'}));
        [td,~] = getModel(td,struct('model_type',model_type,...
            'model_name','S1_joint_smooth','in_signals',...
            {{'opensim',joint_idx}},...
            'out_signals',{'S1_spikes'}));
        [td_act,~] = getModel(td(active_trials),struct('model_type',model_type,...
            'model_name','S1_joint_act','in_signals',...
            {{'opensim',joint_idx}},...
            'out_signals',{'S1_spikes'}));
        [td_pas,~] = getModel(td(passive_trials),struct('model_type',model_type,...
            'model_name','S1_joint_pas','in_signals',...
            {{'opensim',joint_idx}},...
            'out_signals',{'S1_spikes'}));
        % [td,model_info_act] = getModel(td,struct('model_type',model_type,...
        %     'model_name','S1_muscle_act','in_signals',...
        %     {{'opensim_pca',1:6}},...
        %     'out_signals',{'S1_FR'},'train_idx',active_trials));
        % [td,model_info_pas] = getModel(td,struct('model_type',model_type,...
        %     'model_name','S1_muscle_pas','in_signals',...
        %     {{'opensim_pca',1:6}},...
        %     'out_signals',{'S1_FR'},'train_idx',passive_trials));
        % [td_act,model_info_act] = getModel(td(active_trials),struct('model_type',model_type,...
        %     'model_name','S1_muscle_combined','in_signals',...
        %     {{'opensim_pca',1:6}},...
        %     'out_signals',{'S1_FR'}));
        % [td_pas,model_info_pas] = getModel(td(passive_trials),struct('model_type',model_type,...
        %     'model_name','S1_muscle_combined','in_signals',...
        %     {{'opensim_pca',1:6}},...
        %     'out_signals',{'S1_FR'}));
        % td = [td_act,td_pas];
        
        % get dpca for muscle-based neurons
        [~,musc_dpca] = getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{[model_type '_S1_muscle_smooth']}},'num_dims',8,'do_plot',true,...
                                                                            'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                            'dpca_plot_fcn',@dpca_plot_actpas));
        % project into original dPCA space
        getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{[model_type '_S1_muscle_smooth']}},'num_dims',8,'do_plot',true,...
                                                                            'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
                                                                            'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));
        % % get dpca for active-trained muscle-based neurons, projected into original space
        % getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{[model_type 'S1_muscle_act']}},'num_dims',8,'do_plot',true,...
        %                                                                     'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
        %                                                                     'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));
        % % get dpca for passive-trained muscle-based neurons, projected into original space
        % getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{[model_type 'S1_muscle_pas']}},'num_dims',8,'do_plot',true,...
        %                                                                     'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
        %                                                                     'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));
        % % get dpca for combined active-trained and passive-trained muscle-based neurons, projected into original space
        % getDPCA(td,'target_direction','ctrHoldBump',struct('signals',{{[model_type 'S1_muscle_combined']}},'num_dims',8,'do_plot',true,...
        %                                                                     'marg_names',{{'time','direction','active/passive','direction/dynamics interaction'}},...
        %                                                                     'dpca_plot_fcn',@dpca_plot_actpas,'W',dPCA_info.W,'V',dPCA_info.V,'which_marg',dPCA_info.which_marg));

    %% Plot average neural firing rates against muscle-based predictions
        % split td and average
        % td_avg = trialAverage(td,{'ctrHoldBump','target_direction'});
        % [~,td_act] = getTDidx(td_avg,'ctrHoldBump',false);
        % [~,td_pas] = getTDidx(td_avg,'ctrHoldBump',true);
        td_act_avg = trialAverage(td_act,{'ctrHoldBump','target_direction'});
        td_pas_avg = trialAverage(td_pas,{'ctrHoldBump','target_direction'});

        % plot for each neuron
        dir_colors = linspecer(4);
        for neuron_idx = 1:size(td(1).S1_spikes,2)
            h=figure;
            % ylim = [0 max(get_vars(td_avg,{'S1_spikes',neuron_idx}))];
            for dir_idx = 1:length(td_act_avg)
                ax((dir_idx-1)*2+1) = subplot(4,2,(dir_idx-1)*2+1);
                plot(td_act_avg(dir_idx).S1_spikes(:,neuron_idx),'-','linewidth',2,'color',dir_colors(dir_idx,:))
                hold on
                plot(td_act_avg(dir_idx).nn_S1_joint_smooth(:,neuron_idx),'-.','linewidth',2,'color',dir_colors(dir_idx,:))
                plot(td_act_avg(dir_idx).nn_S1_joint_act(:,neuron_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
                % plot(td_act(dir_idx).S1_FR(:,neuron_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
                axis tight
                set(gca,'box','off','tickdir','out')

                ax((dir_idx-1)*2+2) = subplot(4,2,(dir_idx-1)*2+2);
                plot(td_pas_avg(dir_idx).S1_spikes(:,neuron_idx),'-','linewidth',2,'color',dir_colors(dir_idx,:))
                hold on
                plot(td_pas_avg(dir_idx).nn_S1_joint_smooth(:,neuron_idx),'-.','linewidth',2,'color',dir_colors(dir_idx,:))
                plot(td_pas_avg(dir_idx).nn_S1_joint_pas(:,neuron_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
                % plot(td_pas(dir_idx).S1_FR(:,neuron_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
                axis tight
                set(gca,'box','off','tickdir','out')
            end
            linkaxes(ax,'xy')
            waitfor(h)
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Neuron to neuron correlation methods
%% Plot many neurons in active and passive to see whether correlation is broken
    % get a random neuron order
    neuron_order = randperm(size(td(1).S1_spikes,2));

    % split td and average
    td_avg = trialAverage(td,{'ctrHoldBump','target_direction'});
    [~,td_act] = getTDidx(td_avg,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td_avg,'ctrHoldBump',true);

    % maximally distinguishable colors for directions
    dir_colors = linspecer(4);

    % plot each frame with four neurons
    neurons_per_frame = 4;
    for i = 1:neurons_per_frame:length(neuron_order)
        h=figure;
        for j = 1:neurons_per_frame
            % exit loop if we're done with neurons
            if i+j-1>length(neuron_order)
                break;
            end

            % get neuron index
            neuron_idx = neuron_order(i+j-1);
            ylim = [0 max(get_vars(td_avg,{'S1_spikes',neuron_idx}))];

            % active
            subplot(4,2,(j-1)*2+1)
            for dir_idx = 1:length(td_act)
                plot(td_act(dir_idx).S1_spikes(:,neuron_idx),'-','linewidth',2,'color',dir_colors(dir_idx,:))
                hold on
            end
            axis tight
            set(gca,'box','off','tickdir','out','ylim',ylim)
            title(sprintf('Neuron %d - Active',neuron_idx))

            % passive
            subplot(4,2,(j-1)*2+2)
            for dir_idx = 1:length(td_pas)
                plot(td_pas(dir_idx).S1_spikes(:,neuron_idx),'-','linewidth',2,'color',dir_colors(dir_idx,:))
                hold on
            end
            axis tight
            set(gca,'box','off','tickdir','out','ylim',ylim)
            title(sprintf('Neuron %d - Passive',neuron_idx))
        end
        waitfor(h)
    end

%% Compare PC spaces of active and passive movements
    % split td and trim to only post movement/bump epoch
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});

    % get PCA for both only on movement portions
    [~,pca_act] = getPCA(td_act,struct('signals',{{'S1_spikes'}},'do_plot',true));
    [~,pca_pas] = getPCA(td_pas,struct('signals',{{'S1_spikes'}},'do_plot',true));
    
    % split td again to recover whole signal
    % [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    % [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

    % project full signal into pc space
    td_act = getPCA(td_act,pca_act);
    td_pas = getPCA(td_pas,pca_pas);

    % project from active into passive pc space and vice versa
    td_cross_act = getPCA(td_act,pca_pas);
    td_cross_pas = getPCA(td_pas,pca_act);

    % plot first two PCs of both
    % td_act_avg = trialAverage(td_act,'target_direction');
    % dir_colors = linspecer(length(td_act_avg));
    % figure
    % for dir_idx = 1:length(td_act_avg)
    %     plot3(td_act_avg(dir_idx).S1_pca(:,1),td_act_avg(dir_idx).S1_pca(:,2),td_act_avg(dir_idx).S1_pca(:,3),'-','linewidth',2,'color',dir_colors(dir_idx,:))
    %     hold on
    %     plot3(td_act_avg(dir_idx).S1_pca(end,1),td_act_avg(dir_idx).S1_pca(end,2),td_act_avg(dir_idx).S1_pca(end,3),'.','markersize',30,'color',dir_colors(dir_idx,:))
    % end
    % % passive
    % td_pas_avg = trialAverage(td_pas,'target_direction');
    % figure
    % for dir_idx = 1:length(td_pas_avg)
    %     plot3(td_pas_avg(dir_idx).S1_pca(:,1),td_pas_avg(dir_idx).S1_pca(:,2),td_pas_avg(dir_idx).S1_pca(:,3),'-','linewidth',2,'color',dir_colors(dir_idx,:))
    %     hold on
    %     plot3(td_pas_avg(dir_idx).S1_pca(end,1),td_pas_avg(dir_idx).S1_pca(end,2),td_pas_avg(dir_idx).S1_pca(end,3),'.','markersize',30,'color',dir_colors(dir_idx,:))
    % end
    % % plot accross predictions
    % td_cross_act_avg = trialAverage(td_cross_act,'target_direction');
    % figure
    % for dir_idx = 1:length(td_cross_act_avg)
    %     plot3(td_cross_act_avg(dir_idx).S1_pca(:,1),td_cross_act_avg(dir_idx).S1_pca(:,2),td_cross_act_avg(dir_idx).S1_pca(:,3),'-','linewidth',2,'color',dir_colors(dir_idx,:))
    %     hold on
    %     plot3(td_cross_act_avg(dir_idx).S1_pca(end,1),td_cross_act_avg(dir_idx).S1_pca(end,2),td_cross_act_avg(dir_idx).S1_pca(end,3),'.','markersize',30,'color',dir_colors(dir_idx,:))
    % end
    % % passive
    % td_cross_pas_avg = trialAverage(td_cross_pas,'target_direction');
    % figure
    % for dir_idx = 1:length(td_cross_pas_avg)
    %     plot3(td_cross_pas_avg(dir_idx).S1_pca(:,1),td_cross_pas_avg(dir_idx).S1_pca(:,2),td_cross_pas_avg(dir_idx).S1_pca(:,3),'-','linewidth',2,'color',dir_colors(dir_idx,:))
    %     hold on
    %     plot3(td_cross_pas_avg(dir_idx).S1_pca(end,1),td_cross_pas_avg(dir_idx).S1_pca(end,2),td_cross_pas_avg(dir_idx).S1_pca(end,3),'.','markersize',30,'color',dir_colors(dir_idx,:))
    % end

    % Variance explained by PCs plot
    S1_pca_act = get_vars(td_act,{'S1_pca',1:size(td_act(1).S1_pca,2)});
    S1_pca_pas = get_vars(td_pas,{'S1_pca',1:size(td_pas(1).S1_pca,2)});
    S1_pca_cross_act = get_vars(td_cross_act,{'S1_pca',1:size(td_cross_act(1).S1_pca,2)});
    S1_pca_cross_pas = get_vars(td_cross_pas,{'S1_pca',1:size(td_cross_pas(1).S1_pca,2)});

    var_act = var(S1_pca_act);
    var_act = var_act/sum(var_act);
    var_pas = var(S1_pca_pas);
    var_pas = var_pas/sum(var_pas);
    var_cross_act = var(S1_pca_cross_act);
    var_cross_act = var_cross_act/sum(var_cross_act);
    var_cross_pas = var(S1_pca_cross_pas);
    var_cross_pas = var_cross_pas/sum(var_cross_pas);
    align_index_pas = sum(var_cross_pas(1:10))/sum(var_pas(1:10));
    align_index_act = sum(var_cross_act(1:10))/sum(var_act(1:10));

    figure
    subplot(2,1,1)
    bar([var_act;var_cross_act]')
    set(gca,'box','off','tickdir','out')
    ylabel '% VAF'
    title 'Active subspace projection VAF by principle component'
    subplot(2,1,2)
    bar([var_pas;var_cross_pas]')
    set(gca,'box','off','tickdir','out')
    ylabel '% VAF'
    xlabel 'Component number'
    title 'Passive subspace projection VAF by principle component'


    figure
    subplot(2,1,1)
    plot(cumsum([var_act;var_cross_act]'),'linewidth',2)
    set(gca,'box','off','tickdir','out','ylim',[0 1],'xlim',[1 length(var_act)])
    ylabel '% VAF'
    title 'Passive subspace projection cumulative VAF by principle component'
    subplot(2,1,2)
    plot(cumsum([var_pas;var_cross_pas]'),'linewidth',2)
    set(gca,'box','off','tickdir','out','ylim',[0 1],'xlim',[1 length(var_act)])
    ylabel '% VAF'
    xlabel 'Component number'
    title 'Passive subspace projection cumulative VAF by principle component'

    % get principle angles between subspaces
    num_dim = 7;
    basis_act = pca_act.w(:,1:num_dim);
    basis_pas = pca_pas.w(:,1:num_dim);
    
    % Get principle angles with SVD
    neural_angles = principal_angles(basis_act,basis_pas);

    % bootstrap to get angles for within conditions
    num_boots = 1000;
    act_angles_boot = zeros(num_dim,num_boots);
    pas_angles_boot = zeros(num_dim,num_boots);
    for bootnum = 1:num_boots
        % sample trials from td_act and td_pas
        [idx_a,idx_b] = crossvalind('holdout',length(td_act),0.5);
        td_act_a = td_act(idx_a);
        td_act_b = td_act(idx_b);
        [idx_a,idx_b] = crossvalind('holdout',length(td_pas),0.5);
        td_pas_a = td_pas(idx_a);
        td_pas_b = td_pas(idx_b);

        % calculate PCA
        [~,pca_act_a] = getPCA(td_act_a,struct('signals',{{'S1_spikes'}},'do_plot',false));
        [~,pca_act_b] = getPCA(td_act_b,struct('signals',{{'S1_spikes'}},'do_plot',false));
        [~,pca_pas_a] = getPCA(td_pas_a,struct('signals',{{'S1_spikes'}},'do_plot',false));
        [~,pca_pas_b] = getPCA(td_pas_b,struct('signals',{{'S1_spikes'}},'do_plot',false));

        % get basis vectors
        basis_act_a = pca_act_a.w(:,1:num_dim);
        basis_act_b = pca_act_b.w(:,1:num_dim);
        basis_pas_a = pca_pas_a.w(:,1:num_dim);
        basis_pas_b = pca_pas_b.w(:,1:num_dim);

        % calculte principal angles
        act_angles_boot(:,bootnum) = principal_angles(basis_act_a,basis_act_b);
        pas_angles_boot(:,bootnum) = principal_angles(basis_pas_a,basis_pas_b);
    end

    act_angles_CI = prctile(act_angles_boot,[2.5 97.5],2);
    pas_angles_CI = prctile(pas_angles_boot,[2.5 97.5],2);
    act_angles_patch = [act_angles_CI(:,1); flipud(act_angles_CI(:,2))];
    pas_angles_patch = [pas_angles_CI(:,1); flipud(pas_angles_CI(:,2))];
    num_patch = [1:7 7:-1:1]';
    
    % plot angles
    figure
    p = patch(num_patch,act_angles_patch,'b');
    set(p,'facealpha',0.5)
    hold on
    p = patch(num_patch,pas_angles_patch,'y');
    set(p,'facealpha',0.5)
    plot(1:7,neural_angles,'-k','linewidth',2)
    set(gca,'box','off','tickdir','out')
    xlabel 'Component number'
    ylabel 'Angle (rad)'
    title 'Principal angles between active and passive'
    legend('Within Active','Within Passive','Active vs. Passive')

    % Calculate angles between first two PCs of active and passive
    crit_angle = acosd(3.3/sqrt(size(basis_act,1)));
    angle1to1 = acosd(basis_act(:,1)'*basis_pas(:,1));
    angle2to2 = acosd(basis_act(:,2)'*basis_pas(:,2));
    angle1to2 = acosd(basis_act(:,1)'*basis_pas(:,2));
    angle2to1 = acosd(basis_act(:,2)'*basis_pas(:,1));

%% Elsayed-style correlation space plots
    %% get pairwise corr in the two epochs
    % split td and trim to only post movement/bump epoch
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});

    % get correlations
    [rho_act,sort_idx] = pairwiseCorr(td_act,struct('signals',{{'S1_spikes'}},'cluster_order',true));
    [rho_pas] = pairwiseCorr(td_pas,struct('signals',{{'S1_spikes'}}));
    % get correlation of correlation
    S1_corr_corr = corr2(rho_act,rho_pas);
    figure
    subplot(2,2,1)
    imagesc(rho_act(sort_idx,sort_idx))
    clim = get(gca,'clim');
    colorbar
    colormap(cm_viridis);
    axis square
    title 'S1 neural correlation - Active'
    subplot(2,2,3)
    imagesc(rho_pas(sort_idx,sort_idx),clim)
    colorbar
    colormap(cm_viridis);
    axis square
    title 'S1 neural correlation - Passive'
    xlabel(sprintf('Similarity index: %f',S1_corr_corr))

    % plot correlations against each other
    subplot(2,2,2)
    scatter(rho_act(:),rho_pas(:),[],'k','filled')
    set(gca,'box','off','tickdir','out','xlim',[-1 1],'ylim',[-1 1])
    title(sprintf('Active-passive corrlation relationship: R^2 = %f',S1_corr_corr^2))
    xlabel 'S1 neural correlation - Active'
    ylabel 'S1 neural correlation - Passive'
    axis equal

    % get epoch preference index
    num_neurons = size(td(1).S1_spikes,2);
    S1_act = get_vars(td_act,{'S1_spikes',1:num_neurons});
    S1_pas = get_vars(td_pas,{'S1_spikes',1:num_neurons});
    % calculate tuning to each epoch
    act_tuning = range(S1_act,1)./mean([S1_act;S1_pas],1);
    pas_tuning = range(S1_pas,1)./mean([S1_pas;S1_pas],1);
    % calculate index by normalizing tuning by average tuning in each epoch and subtract
    pref_index = pas_tuning/mean(pas_tuning) - act_tuning/mean(act_tuning);
    % plot histogram
    subplot(2,2,4)
    histogram(pref_index)
    set(gca,'box','off','tickdir','out')

    %% correlation analysis in muscle space
    % get covariates to base simulated neurons on
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));

    % get active and passive trial indices
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

    % correlations and meta-correlations
    [rho_act,sort_idx] = pairwiseCorr(td_act,struct('signals',{{'opensim',opensim_idx}},'cluster_order',true));
    [rho_pas] = pairwiseCorr(td_pas,struct('signals',{{'opensim',opensim_idx}}));
    % get correlation of correlation
    musc_corr_corr = corr2(rho_act,rho_pas);
    figure
    subplot(2,2,1)
    imagesc(rho_act(sort_idx,sort_idx))
    clim = get(gca,'clim');
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Muscle velocity correlation - Active'
    subplot(2,2,3)
    imagesc(rho_pas(sort_idx,sort_idx),clim)
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Muscle velocity correlation - Passive'
    xlabel(sprintf('Similarity index: %f',musc_corr_corr))

    % plot correlations against each other
    subplot(2,2,2)
    scatter(rho_act(:),rho_pas(:),[],'k','filled')
    set(gca,'box','off','tickdir','out','xlim',[-1 1],'ylim',[-1 1])
    title(sprintf('Active-passive corrlation relationship: R^2 = %f',musc_corr_corr^2))
    xlabel 'Muscle neural correlation - Active'
    ylabel 'Muscle neural correlation - Passive'
    axis equal

    % get epoch preference index
    musc_act = get_vars(td_act,{'opensim',opensim_idx});
    musc_pas = get_vars(td_pas,{'opensim',opensim_idx});
    % calculate tuning to each epoch
    act_tuning = range(musc_act,1)./mean([musc_act;musc_pas],1);
    pas_tuning = range(musc_pas,1)./mean([musc_pas;musc_pas],1);
    % calculate index by normalizing tuning by average tuning in each epoch and subtract
    pref_index = pas_tuning/mean(pas_tuning) - act_tuning/mean(act_tuning);
    % plot histogram
    subplot(2,2,4)
    histogram(pref_index)
    set(gca,'box','off','tickdir','out')

    %% correlation analysis in joint space
    % get covariates to base simulated neurons on
    opensim_idx = find(contains(td(1).opensim_names,'_vel'));

    % get active and passive trial indices
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

    % correlations and meta correlations
    [rho_act,sort_idx] = pairwiseCorr(td_act,struct('signals',{{'opensim',opensim_idx}},'cluster_order',true));
    [rho_pas] = pairwiseCorr(td_pas,struct('signals',{{'opensim',opensim_idx}}));
    % get correlation of correlation
    joint_corr_corr = corr2(rho_act,rho_pas);
    figure
    subplot(2,2,1)
    imagesc(rho_act(sort_idx,sort_idx))
    clim = get(gca,'clim');
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Joint velocity correlation - Active'
    subplot(2,2,3)
    imagesc(rho_pas(sort_idx,sort_idx),clim)
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Joint velocity correlation - Passive'
    xlabel(sprintf('Similarity index: %f',joint_corr_corr))

    % plot correlations against each other
    subplot(2,2,2)
    scatter(rho_act(:),rho_pas(:),[],'k','filled')
    set(gca,'box','off','tickdir','out','xlim',[-1 1],'ylim',[-1 1])
    title(sprintf('Active-passive corrlation relationship: R^2 = %f',joint_corr_corr^2))
    xlabel 'Joint neural correlation - Active'
    ylabel 'Joint neural correlation - Passive'
    axis equal

    % get epoch preference index
    joint_act = get_vars(td_act,{'opensim',opensim_idx});
    joint_pas = get_vars(td_pas,{'opensim',opensim_idx});
    % calculate tuning to each epoch
    act_tuning = range(joint_act,1)./mean([joint_act;joint_pas],1);
    pas_tuning = range(joint_pas,1)./mean([joint_pas;joint_pas],1);
    % calculate index by normalizing tuning by average tuning in each epoch and subtract
    pref_index = pas_tuning/mean(pas_tuning) - act_tuning/mean(act_tuning);
    % plot histogram
    subplot(2,2,4)
    histogram(pref_index)
    set(gca,'box','off','tickdir','out')

    %% correlation analysis of muscle based neurons
    % get covariates to base simulated neurons on
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
    td = getPCA(td,struct('signals',{{'opensim',opensim_idx}},'do_plot',false));

    % split td and trim to only post movement/bump epoch
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});

    % correlation and meta correlation
    [rho_act,sort_idx] = pairwiseCorr(td_act,struct('signals',{{[model_type '_S1_muscle']}},'cluster_order',true));
    [rho_pas] = pairwiseCorr(td_pas,struct('signals',{{[model_type '_S1_muscle']}}));
    % meta correlation
    musc_model_corr_corr = corr2(rho_act,rho_pas);
    % plot
    figure
    subplot(2,2,1)
    imagesc(rho_act(sort_idx,sort_idx))
    clim = get(gca,'clim');
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Muscle-based neural correlation - Active'
    subplot(2,2,3)
    imagesc(rho_pas(sort_idx,sort_idx),clim)
    colorbar
    colormap(cm_viridis);
    axis square
    title 'Muscle-based neural correlation - Passive'
    xlabel(sprintf('Similarity index: %f',musc_model_corr_corr))

    % plot correlations against each other
    subplot(2,2,2)
    scatter(rho_act(:),rho_pas(:),[],'k','filled')
    set(gca,'box','off','tickdir','out','xlim',[-1 1],'ylim',[-1 1])
    title(sprintf('Active-passive corrlation relationship: R^2 = %f',musc_model_corr_corr^2))
    xlabel 'Muscle-based neural correlation - Active'
    ylabel 'Muscle-based neural correlation - Passive'

    % get epoch preference index
    num_neurons = size(td(1).S1_spikes,2);
    S1_act = get_vars(td_act,{[model_type '_S1_muscle'],1:num_neurons});
    S1_pas = get_vars(td_pas,{[model_type '_S1_muscle'],1:num_neurons});
    % calculate tuning to each epoch
    act_tuning = range(S1_act,1)./mean([S1_act;S1_pas],1);
    pas_tuning = range(S1_pas,1)./mean([S1_pas;S1_pas],1);
    % calculate index by normalizing tuning by average tuning in each epoch and subtract
    pref_index = pas_tuning/mean(pas_tuning) - act_tuning/mean(act_tuning);
    % plot histogram
    subplot(2,2,4)
    histogram(pref_index)
    set(gca,'box','off','tickdir','out')
    
    %% Test pairwise correlation stability in each epoch
    % split td and trim to only post movement/bump epoch
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});

    % which one to test...
    td_test = td_pas;
    test_var = 'Passive';

    % crossval split
    [idx_1,idx_2] = crossvalind('HoldOut',length(td_test),0.5);

    [rho_1,sort_idx] = pairwiseCorr(td_test(idx_1),struct('signals',{{'S1_spikes'}},'cluster_order',true));
    [rho_2] = pairwiseCorr(td_test(idx_2),struct('signals',{{'S1_spikes'}}));

    S1_corr_corr = corr2(rho_1,rho_2);

    figure
    subplot(2,2,1)
    imagesc(rho_1(sort_idx,sort_idx))
    clim = get(gca,'clim');
    colorbar
    colormap(cm_viridis);
    axis square
    title(sprintf('S1 neural correlation - %s',test_var))
    subplot(2,2,3)
    imagesc(rho_2(sort_idx,sort_idx),clim)
    colorbar
    colormap(cm_viridis);
    axis square
    title(sprintf('S1 neural correlation - %s',test_var))
    xlabel(sprintf('Similarity index: %f',S1_corr_corr))

    % plot correlations against each other
    subplot(2,2,2)
    scatter(rho_1(:),rho_2(:),[],'k','filled')
    set(gca,'box','off','tickdir','out','xlim',[-1 1],'ylim',[-1 1])
    title(sprintf('%s-%s corrlation relationship: R^2 = %f',test_var,test_var,S1_corr_corr^2))
    xlabel(sprintf('S1 neural correlation - %s',test_var))
    ylabel(sprintf('S1 neural correlation - %s',test_var))

%% Look at canonical correlations between latent neural space and behavior in active and passive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basics...
    %% check spike means and variability between act and pas
        % split td and trim to only post movement/bump epoch
        [~,td_act] = getTDidx(td,'ctrHoldBump',false);
        td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
        [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
        td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});

        spikes_act = get_vars(td_act,check_signals(td_act,'S1_spikes'));
        spikes_pas = get_vars(td_pas,check_signals(td_pas,'S1_spikes'));

        means_act = mean(spikes_act);
        means_pas = mean(spikes_pas);
        var_act = var(spikes_act);
        var_pas = var(spikes_pas);

        figure
        subplot(2,1,1)
        bar([means_act;means_pas]')
        set(gca,'box','off','tickdir','out')
        title 'Firing rate mean'
        ylabel 'Hz'
        xlabel 'Neuron Number'
        legend('Active','Passive')
        subplot(2,1,2)
        bar([var_act;var_pas]')
        set(gca,'box','off','tickdir','out')
        title 'Firing rate variance'
        ylabel 'Hz^2'
        xlabel 'Neuron Number'
        legend('Active','Passive')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Misc Visualizations
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

%% Plot velocity traces
    % split td
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

    % get directions
    dirs = sort(unique([td.target_direction]));
    dir_colors = linspecer(length(dirs));

    % get mean and CI of velocity traces (saved in tensors)
    mean_vel_act = zeros(size(td(1).vel,1),size(td(1).vel,2),length(dirs));
    CI_vel_act = zeros(size(td(1).vel,1),size(td(1).vel,2),length(dirs),2);
    mean_vel_pas = zeros(size(td(1).vel,1),size(td(1).vel,2),length(dirs));
    CI_vel_pas = zeros(size(td(1).vel,1),size(td(1).vel,2),length(dirs),2);
    for dir_idx = 1:length(dirs)
        % active
        [~,td_temp] = getTDidx(td_act,'target_direction',dirs(dir_idx));
        vel_act = cat(3,td_temp.vel);
        mean_vel_act(:,:,dir_idx) = mean(vel_act,3);
        stderr = std(vel_act,0,3)/sqrt(size(vel_act,3));
        t_lo = tinv(0.025,size(vel_act,3)-1);
        t_hi = tinv(0.975,size(vel_act,3)-1);
        CI_vel_act(:,:,dir_idx,1) = mean_vel_act(:,:,dir_idx)-stderr*t_lo;
        CI_vel_act(:,:,dir_idx,2) = mean_vel_act(:,:,dir_idx)-stderr*t_hi;

        % passive
        [~,td_temp] = getTDidx(td_pas,'target_direction',dirs(dir_idx));
        vel_pas = cat(3,td_temp.vel);
        mean_vel_pas(:,:,dir_idx) = mean(vel_pas,3);
        stderr = std(vel_pas,0,3)/sqrt(size(vel_pas,3));
        t_lo = tinv(0.025,size(vel_pas,3)-1);
        t_hi = tinv(0.975,size(vel_pas,3)-1);
        CI_vel_pas(:,:,dir_idx,1) = mean_vel_pas(:,:,dir_idx)-stderr*t_lo;
        CI_vel_pas(:,:,dir_idx,2) = mean_vel_pas(:,:,dir_idx)-stderr*t_hi;
    end

    % plot it out
    figure
    % ax = zeros(4*2);
    ax = zeros(2,1);
    for dir_idx = 1:length(dirs)
        for sig_idx = 1:2
            % plt_idx = (dir_idx-1)*2+sig_idx;
            % ax(plt_idx) = subplot(4,2,plt_idx);
            ax(sig_idx) = subplot(2,1,sig_idx);

            % active
            patchtrace_x = [1:size(CI_vel_act,1) size(CI_vel_act,1):-1:1]';
            patchtrace_y = cat(1,CI_vel_act(:,sig_idx,dir_idx,1),CI_vel_act(end:-1:1,sig_idx,dir_idx,2));
            patch(patchtrace_x,patchtrace_y,dir_colors(dir_idx,:),'facealpha',0.5,'edgecolor','none')
            hold on
            plot(mean_vel_act(:,sig_idx,dir_idx),'-','linewidth',2,'color',dir_colors(dir_idx,:))

            % passive
            patchtrace_x = [1:size(CI_vel_pas,1) size(CI_vel_pas,1):-1:1]';
            patchtrace_y = cat(1,CI_vel_pas(:,sig_idx,dir_idx,1),CI_vel_pas(end:-1:1,sig_idx,dir_idx,2));
            patch(patchtrace_x,patchtrace_y,dir_colors(dir_idx,:),'facealpha',0.5,'edgecolor','none')
            hold on
            plot(mean_vel_pas(:,sig_idx,dir_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
        end
    end
    linkaxes(ax)

    figure
    timevec = (1:size(mean_vel_act,1))';
    for dir_idx = 1:length(dirs)
        plot3(timevec,mean_vel_act(:,1,dir_idx),mean_vel_act(:,2,dir_idx),'-','linewidth',2,'color',dir_colors(dir_idx,:))
        hold on
        plot3(timevec,mean_vel_pas(:,1,dir_idx),mean_vel_pas(:,2,dir_idx),'--','linewidth',2,'color',dir_colors(dir_idx,:))
    end

    % % plot traces
    % figure
    % for dir_idx = 1:length(dirs)
    %     [~,td_temp] = getTDidx(td_pas,'target_direction',dirs(dir_idx));
    %     for i = 1:length(td_temp)
    %         subplot(4,2,(dir_idx-1)*2+1)
    %         plot(td_temp(i).vel(:,1),'-','linewidth',1,'color',dir_colors(dir_idx,:))
    %         hold on

    %         subplot(4,2,(dir_idx-1)*2+2)
    %         plot(td_temp(i).vel(:,2),'-','linewidth',1,'color',dir_colors(dir_idx,:))
    %         hold on
    %     end
    % end

    % % plot traces
    % figure
    % for dir_idx = 1:length(dirs)
    %     [~,td_temp] = getTDidx(td_act,'target_direction',dirs(dir_idx));
    %     td_temp = trimTD(td_temp,{'idx_movement_on',0},{'idx_movement_on',12});
    %     subplot(2,1,1)
    %     for i = 1:length(td_temp)
    %         plot(td_temp(i).vel(:,1),td_temp(i).vel(:,2),'-','linewidth',1,'color',dir_colors(dir_idx,:))
    %         hold on
    %     end
    %     axis square
    %     axis equal

    %     [~,td_temp] = getTDidx(td_pas,'target_direction',dirs(dir_idx));
    %     td_temp = trimTD(td_temp,{'idx_bumpTime',0},{'idx_bumpTime',12});
    %     subplot(2,1,2)
    %     for i = 1:length(td_temp)
    %         plot(td_temp(i).vel(:,1),td_temp(i).vel(:,2),'-','linewidth',1,'color',dir_colors(dir_idx,:))
    %         hold on
    %     end
    %     axis square
    %     axis equal
    % end

