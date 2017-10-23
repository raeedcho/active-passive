%% set up TD
clearvars -except trial_data
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);

[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});

td = cat(2,td_act,td_pas);

% get pairwise correlations of neurons
[rho_act,sort_idx_act] = pairwiseCorr(td_act,struct('signals',{{'S1_spikes'}},'cluster_order',true));
[rho_pas,sort_idx_pas] = pairwiseCorr(td_pas,struct('signals',{{'S1_spikes'}}));

% plot neuron pairwise correlations for active and passive
clim = [min(min([rho_act rho_pas])) max(max([rho_act rho_pas]))];
figure
subplot(1,2,1)
imagesc(rho_act,clim)
subplot(1,2,2)
imagesc(rho_pas(sort_idx_act,sort_idx_act),clim);

%% try looking at canonical correlations of PC space
% get PCA representations
[~,pca_info_act] = getPCA(td_act,struct('signals',{{'S1_spikes'}},'do_plot',true));
[~,pca_info_pas] = getPCA(td_pas,struct('signals',{{'S1_spikes'}},'do_plot',true));

% extract basis sets
num_dim = 7;
basis_act = pca_info_act.w(:,1:num_dim);
basis_pas = pca_info_pas.w(:,1:num_dim);

% Get canonical correlations with SVD
[~,s,~] = svd(basis_act'*basis_pas);
neural_angles = acosd(diag(s))

%% try looking at canonical correlations of muscle space
% get PCA representations
opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
[~,pca_info_act] = getPCA(td_act,struct('signals',{{'opensim',opensim_idx}},'do_plot',true));
[~,pca_info_pas] = getPCA(td_pas,struct('signals',{{'opensim',opensim_idx}},'do_plot',true));

% extract basis sets
num_dim = 7;
basis_act = pca_info_act.w(:,1:num_dim);
basis_pas = pca_info_pas.w(:,1:num_dim);

% Get canonical correlations with SVD
[~,s,~] = svd(basis_act'*basis_pas);
muscle_angles = acosd(diag(s))

%% Look at potent space from neurons to muscle
[~,pca_info_act] = getPotentSpace(td_act,struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',6,'out_dims',3));
[~,pca_info_pas] = getPotentSpace(td_pas,struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',6,'out_dims',3));

% extract basis sets
basis_act = pca_info_act.V_potent;
basis_pas = pca_info_pas.V_potent;

% Get canonical correlations with SVD
[~,s,~] = svd(basis_act'*basis_pas);
potent_angles = acosd(diag(s));

%% check to see if angle to same subspace is large
opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
[idx1,idx2] = crossvalind('HoldOut',length(td_act),0.5);
[~,pca_info1] = getPotentSpace(td_act(idx1),struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',6,'out_dims',3));
[~,pca_info2] = getPotentSpace(td_act(idx2),struct('in_signals',{{'S1_spikes'}},'out_signals',{{'opensim',opensim_idx}},'in_dims',6,'out_dims',3));

% extract basis sets
basis1 = pca_info1.V_potent;
basis2 = pca_info2.V_potent;

% Get canonical correlations with SVD
[~,s,~] = svd(basis1'*basis2);
self_angles = acosd(diag(s));

%% bootstrap potent space angles
num_boot = 1000
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
