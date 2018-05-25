%% set up TD
    [~,td] = getTDidx(trial_data,'result','R');
    
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
    td = removeBadNeurons(td,struct('do_shunt_check',true));

    td = smoothSignals(td,struct('signals','S1_spikes'));
    
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    td_act = trimTD(td_act,{'idx_movement_on',-10},{'idx_movement_on',15});
    nanners = isnan(cat(1,td_act.target_direction));
    td_act = td_act(~nanners);
    % td_act = binTD(td_act,5);
    
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',-10},{'idx_bumpTime',15});
    % td_pas = binTD(td_pas,5);
    
    td = cat(2,td_act,td_pas);
    
%% Check stats of PCA
    % [td_pas,pca_info_pas] = getPCA(td_pas,struct('signals',{{'S1_spikes'}},'do_plot',true));
    signals = get_vars(td_pas,{'S1_pca',1:12});

%% get pairwise correlations of neurons
    [rho_act,sort_idx_act] = pairwiseCorr(td_act,struct('signals',{{'S1_spikes'}},'cluster_order',true));
    [rho_pas,sort_idx_pas] = pairwiseCorr(td_pas,struct('signals',{{'S1_spikes'}}));
    
    % plot neuron pairwise correlations for active and passive
    clim = [min(min([rho_act rho_pas])) max(max([rho_act rho_pas]))];
    figure
    subplot(1,2,1)
    imagesc(rho_act,clim)
    axis square
    subplot(1,2,2)
    imagesc(rho_pas(sort_idx_act,sort_idx_act),clim);
    axis square

%% try looking at principal angles of PC space
    % get PCA representations
    [~,pca_info_act] = getPCA(td_act,struct('signals',{{'S1_spikes'}},'do_plot',true));
    [~,pca_info_pas] = getPCA(td_pas,struct('signals',{{'S1_spikes'}},'do_plot',true));
    
    % extract basis sets
    num_dim = 7;
    basis_act = pca_info_act.w(:,1:num_dim);
    basis_pas = pca_info_pas.w(:,1:num_dim);
    
    % Get principle angles with SVD
    [~,s,~] = svd(basis_act'*basis_pas);
    neural_angles = acosd(diag(s));
    
    %% try looking at principle angles of muscle space
    % get PCA representations
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
    [~,pca_info_act] = getPCA(td_act,struct('signals',{{'opensim',opensim_idx}},'do_plot',true));
    [~,pca_info_pas] = getPCA(td_pas,struct('signals',{{'opensim',opensim_idx}},'do_plot',true));
    
    % extract basis sets
    num_dim = 7;
    basis_act = pca_info_act.w(:,1:num_dim);
    basis_pas = pca_info_pas.w(:,1:num_dim);
    
    % Get principal angles with SVD
    [~,s,~] = svd(basis_act'*basis_pas);
    muscle_angles = acosd(diag(s));
    
    %% Look at potent space from neurons to muscle
    [~,pca_info_act] = getPotentSpace(td_act,struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',12,'out_dims',6));
    [~,pca_info_pas] = getPotentSpace(td_pas,struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',12,'out_dims',6));
    
    % extract basis sets
    basis_act = pca_info_act.V_potent;
    basis_pas = pca_info_pas.V_potent;
    
    % Get canonical correlations with SVD
    [~,s,~] = svd(basis_act'*basis_pas);
    potent_angles = acosd(diag(s));
    disp(potent_angles)

%% check to see if angle to same subspace is large
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
    [idx1,idx2] = crossvalind('HoldOut',length(td_act),0.5);
    [~,pca_info1] = getPotentSpace(td_act(idx1),struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',12,'out_dims',6));
    [~,pca_info2] = getPotentSpace(td_act(idx2),struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',12,'out_dims',6));
    
    % extract basis sets
    basis1 = pca_info1.V_potent;
    basis2 = pca_info2.V_potent;
    
    % Get canonical correlations with SVD
    [~,s,~] = svd(basis1'*basis2);
    self_angles = acosd(diag(s));
    disp(self_angles)

%% bootstrap potent space angles
    num_boot = 1000;
    potent_angles = zeros(3,num_boot);
    self_angles = zeros(3,num_boot);
    num_samp = min(length(td_act),length(td_pas));
    warn = warning('query','last');
    warning('off',warn.identifier);
    for i=1:1000
        %% Look at potent space from neurons to muscle
        idx = randi(num_samp,1,num_samp);
        [~,pca_info_act] = getPotentSpace(td_act(idx),struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',6,'out_dims',3));
        [~,pca_info_pas] = getPotentSpace(td_pas(idx),struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',6,'out_dims',3));
    
        % extract basis sets
        basis_act = pca_info_act.V_potent;
        basis_pas = pca_info_pas.V_potent;
    
        % Get canonical correlations with SVD
        [~,s,~] = svd(basis_act'*basis_pas);
        potent_angles(:,i) = acosd(diag(s));
    
        %% check to see if angle to same subspace is large
        opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
        idx1 = randi(num_samp,1,num_samp);
        idx2 = randi(num_samp,1,num_samp);
        [~,pca_info1] = getPotentSpace(td_act(idx1),struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',6,'out_dims',3));
        [~,pca_info2] = getPotentSpace(td_act(idx2),struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',6,'out_dims',3));
    
        % extract basis sets
        basis1 = pca_info1.V_potent;
        basis2 = pca_info2.V_potent;
    
        % Get canonical correlations with SVD
        [~,s,~] = svd(basis1'*basis2);
        self_angles(:,i) = acosd(diag(s));
        
        disp(['Finished iteration ' num2str(i) ' of ' num2str(num_boot)])
    end
    warning('on',warn.identifier);
    
    prctile(potent_angles,[2.5 97.5],2)
    prctile(self_angles,[2.5 97.5],2)

%% Figure out what's in the null space of S1->muscle velocity
    in_dims = 12;
    out_dims = 6;
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
    [td_pas,pca_info_pas] = getPotentSpace(td_pas,struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',in_dims,'out_dims',out_dims));

    pca_info = struct('signals',{{'S1_spikes'}},'w',pca_info_pas.w_in,'mu',pca_info_pas.mu_in);
    td_act_proj = getPCA(td_act,pca_info);
    % add potent and null projections to active TD
    for trial = 1:length(td_act_proj)
        data = td_act_proj(trial).S1_pca;
        data = data(:,1:in_dims);
        td_act_proj(trial).S1opensim_potent = data * pca_info_pas.V_potent;
        td_act_proj(trial).S1opensim_null = data * pca_info_pas.V_null;
    end

    % plot active/passive potent and active/passive null
    td_pas_avg = trialAverage(td_pas,'bumpDir');
    % td_pas_avg = td_pas;
    td_act_avg = trialAverage(td_act_proj,'target_direction');
    dir_colors = linspecer(4);
    for trial = 1:length(td_pas_avg)
        %passives
        dir_idx = td_pas_avg(trial).bumpDir/90 + 1;

        figure(1)
        plot3(td_pas_avg(trial).S1opensim_potent(:,1),td_pas_avg(trial).S1opensim_potent(:,2),td_pas_avg(trial).S1opensim_potent(:,3),'--','linewidth',2,'color',dir_colors(dir_idx,:));
        hold on
        plot3(td_pas_avg(trial).S1opensim_potent(end,1),td_pas_avg(trial).S1opensim_potent(end,2),td_pas_avg(trial).S1opensim_potent(end,3),'o','markersize',25,'color',dir_colors(dir_idx,:));
        
        figure(2)
        plot3(td_pas_avg(trial).S1opensim_null(:,1),td_pas_avg(trial).S1opensim_null(:,2),td_pas_avg(trial).S1opensim_null(:,3),'--','linewidth',2,'color',dir_colors(dir_idx,:));
        hold on
        plot3(td_pas_avg(trial).S1opensim_null(end,1),td_pas_avg(trial).S1opensim_null(end,2),td_pas_avg(trial).S1opensim_null(end,3),'o','markersize',25,'color',dir_colors(dir_idx,:));

        figure(3)
        plot3(td_pas_avg(trial).S1_pca(:,1),td_pas_avg(trial).S1_pca(:,2),td_pas_avg(trial).S1_pca(:,3),'--','linewidth',2,'color',dir_colors(dir_idx,:));
        hold on
        plot3(td_pas_avg(trial).S1_pca(end,1),td_pas_avg(trial).S1_pca(end,2),td_pas_avg(trial).S1_pca(end,3),'o','markersize',25,'color',dir_colors(dir_idx,:));

        %actives
        dir_idx = td_act_avg(trial).target_direction/(pi/2) + 1;
        dir_idx = round(dir_idx);

        figure(1)
        plot3(td_act_avg(trial).S1opensim_potent(:,1),td_act_avg(trial).S1opensim_potent(:,2),td_act_avg(trial).S1opensim_potent(:,3),'linewidth',2,'color',dir_colors(dir_idx,:));
        hold on
        plot3(td_act_avg(trial).S1opensim_potent(end,1),td_act_avg(trial).S1opensim_potent(end,2),td_act_avg(trial).S1opensim_potent(end,3),'.','markersize',50,'color',dir_colors(dir_idx,:));
        
        figure(2)
        plot3(td_act_avg(trial).S1opensim_null(:,1),td_act_avg(trial).S1opensim_null(:,2),td_act_avg(trial).S1opensim_null(:,3),'linewidth',2,'color',dir_colors(dir_idx,:));
        hold on
        plot3(td_act_avg(trial).S1opensim_null(end,1),td_act_avg(trial).S1opensim_null(end,2),td_act_avg(trial).S1opensim_null(end,3),'.','markersize',50,'color',dir_colors(dir_idx,:));

        figure(3)
        plot3(td_act_avg(trial).S1_pca(:,1),td_act_avg(trial).S1_pca(:,2),td_act_avg(trial).S1_pca(:,3),'linewidth',2,'color',dir_colors(dir_idx,:));
        hold on
        plot3(td_act_avg(trial).S1_pca(end,1),td_act_avg(trial).S1_pca(end,2),td_act_avg(trial).S1_pca(end,3),'.','markersize',50,'color',dir_colors(dir_idx,:));
    end

%% What's happening in cokernel (left null space) of muscle velocity -> S1 activity?
    S1_dims = 12;
    musc_dims = 6;
    opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));

    % get neural and muscle latent space
    td_pas = getPCA(td_pas,struct('signals',{{'opensim',opensim_idx}},'do_plot',true));
    [td_pas,pca_info] = getPCA(td_pas,struct('signals',{{'S1_spikes'}},'do_plot',true));
    td_act = getPCA(td_act,pca_info);

    % get projection from muscle latent to neural latent
    model_params = struct('model_type','linmodel','model_name','musc_model',...
                            'in_signals',{{'opensim_pca',1:musc_dims}},'out_signals',{{'S1_pca',1:S1_dims}});
    [td_pas,model_info] = getModel(td_pas,model_params);
    % get projection matrix (S1_dims x musc_dims, such that neural = proj_mat * muscle)
    proj_mat = model_info.b(2:end,:)'; % drop offset weight because input and output data is zero-centered

    % decompose proj_mat for input and output spaces
    [u,s,v] = svd(proj_mat);

    % first musc_dims columns of u are basis for image of musc->S1 transform
    image_basis = u(:,1:musc_dims);
    cokernel_basis = u(:,(musc_dims+1):end);

    % project S1_pca into image and cokernel bases
    for trial = 1:length(td_pas)
        data = td_pas(trial).S1_pca(:,1:S1_dims);
        image = data*image_basis;
        cokernel = data*cokernel_basis;
        td_pas(trial).muscS1_image = image;
        td_pas(trial).muscS1_cokernel = cokernel;
    end

    for trial = 1:length(td_act)
        data = td_act(trial).S1_pca(:,1:S1_dims);
        image = data*image_basis;
        cokernel = data*cokernel_basis;
        td_act(trial).muscS1_image = image;
        td_act(trial).muscS1_cokernel = cokernel;
    end

    % Plot magnitude of image against magnitude of cokernel
    % td_pas_avg = trialAverage(td_pas,'bumpDir');
    td_pas_avg = td_pas;
    figure
    dir_colors = linspecer(4);
    for trial = 1:length(td_pas_avg)
        image = td_pas_avg(trial).muscS1_image;
        image_mag = sqrt(sum(image.^2,2));

        cokernel = td_pas_avg(trial).muscS1_cokernel;
        cokernel_mag = sqrt(sum(cokernel.^2,2));

        % plot it out
        dir_idx = td_pas_avg(trial).bumpDir/90 + 1;
        plot(cokernel_mag,image_mag,'--','linewidth',2,'color',dir_colors(dir_idx,:))
        hold on
        plot(cokernel_mag(end),image_mag(end),'o','markersize',25,'color',dir_colors(dir_idx,:))
    end
    legend({'0','','90','','180','','270'})
    axis equal
    set(gca,'box','off','tickdir','out')

    % Plot magnitude of image against magnitude of cokernel
    % td_act_avg = trialAverage(td_act,'target_direction');
    td_act_avg = td_act;
    %figure
    for trial = 1:length(td_act_avg)
        image = td_act_avg(trial).muscS1_image;
        image_mag = sqrt(sum(image.^2,2));

        cokernel = td_act_avg(trial).muscS1_cokernel;
        cokernel_mag = sqrt(sum(cokernel.^2,2));

        % plot it out
        dir_idx = td_act_avg(trial).target_direction/(pi/2) + 1;
        dir_idx = round(dir_idx);
        plot(cokernel_mag,image_mag,'-','linewidth',2,'color',dir_colors(dir_idx,:))
        hold on
        plot(cokernel_mag(end),image_mag(end),'.','markersize',50,'color',dir_colors(dir_idx,:))
    end
    legend({'0','','90','','180','','270'})
    axis equal
    set(gca,'box','off','tickdir','out')
