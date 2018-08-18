%% Get act vs pas
    [~,td] = getTDidx(trial_data,'result','R');
    
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
    td = smoothSignals(td,struct('signals','markers'));
    td = getDifferential(td,struct('signal','markers','alias','marker_vel'));
    % add firing rates rather than spike counts
    td = addFiringRates(td,struct('array','S1'));
    
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    % clean nans out...?
    nanners = isnan(cat(1,td_act.target_direction));
    td_act = td_act(~nanners);
    td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',14});
    td_act = binTD(td_act,'average');
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',14});
    td_pas = binTD(td_pas,'average');
    td = cat(2,td_act,td_pas);

%% get splits on array
    if exist('ant_chans','var')
        % array_signals = find(ismember(td(1).S1_unit_guide(:,1),ant_chans));
        array_signals = find(~ismember(td(1).S1_unit_guide(:,1),ant_chans));
        % array_signals = sort(randperm(length(td(1).S1_unit_guide(:,1)),floor(length(td(1).S1_unit_guide(:,1))/2))');
        % array_signals = 1:length(td(1).S1_unit_guide);
    else
        array_signals = 1:length(td(1).S1_unit_guide);
    end

%% Get models on data
    hand_idx = 1:3;
    elbow_idx = 28:30;
    % Hand kinematics model
    [td,model_info] = getModel(td,struct('model_type','glm',...
        'model_name','ext','in_signals',{{'markers',hand_idx;'marker_vel',hand_idx}},...
        'out_signals',{{'S1_FR',array_signals}}));

    % Hand forcekin model
    [td,model_info] = getModel(td,struct('model_type','glm',...
        'model_name','forcekin','in_signals',{{'markers',hand_idx;'marker_vel',hand_idx;'force','all'}},...
        'out_signals',{{'S1_FR',array_signals}}));

    % Hand/Elbow model
    [td,model_info] = getModel(td,struct('model_type','glm',...
        'model_name','hand_elbow','in_signals',{{'markers',[hand_idx elbow_idx];'marker_vel',[hand_idx elbow_idx]}},...
        'out_signals',{{'S1_FR',array_signals}}));

%% get PCA
    td = sqrtTransform(td,'S1_FR');
    % td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
    td = dimReduce(td,struct('signals',{{'S1_FR',array_signals}}));
    
    td = sqrtTransform(td,'glm_ext');
    % td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
    td = dimReduce(td,struct('signals',{{'glm_ext',array_signals}}));

    td = sqrtTransform(td,'glm_forcekin');
    % td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
    td = dimReduce(td,struct('signals',{{'glm_forcekin',array_signals}}));

    td = sqrtTransform(td,'glm_hand_elbow');
    % td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
    td = dimReduce(td,struct('signals',{{'glm_hand_elbow',array_signals}}));
    
%% test separability
    figure(1)
    [actual_sep,actual_mdl] = test_sep(td,struct('signals',{{'S1_FR_pca'}},'do_plot',true));
    title(['Actual separability - ' num2str(actual_sep)])
    
    % % then for directional separability/other view
    % [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    % [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    % signal_act = cat(1,td_act.S1_FR_pca);
    % signal_pas = cat(1,td_pas.S1_FR_pca);
    % % plot active as filled, passive as open
    % bump_colors = linspecer(4);
    % act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
    % pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;
    % figure(2)
    % hold all
    % scatter3(signal_act(:,1),signal_act(:,2),signal_act(:,3),50,bump_colors(act_dir_idx,:),'filled')
    % scatter3(signal_pas(:,1),signal_pas(:,2),signal_pas(:,3),100,bump_colors(pas_dir_idx,:),'o','linewidth',2)
    % ylim = get(gca,'ylim');
    % zlim = get(gca,'zlim');
    % set(gca,'box','off','tickdir','out')
    % axis equal
    % axis off
    % title 'directional sep'

%% Try fabricating trial_data with linear models based on handle kinematics and force
    figure(3)
    [ext_sep,ext_mdl] = test_sep(td,struct('signals',{{'glm_ext_pca',1:4}},'do_plot',true));
    title(['ext separability - ' num2str(ext_sep)])
    figure(4)
    [ext_sep,ext_mdl] = test_sep(td,struct('signals',{{'glm_ext_pca'}},'mdl',actual_mdl,'do_plot',true));
    title(['ext w/ actual disc. separability - ' num2str(ext_sep)])

    figure(5)
    [velforce_sep,velforce_mdl] = test_sep(td,struct('signals',{{'glm_forcekin_pca',1:4}},'do_plot',true));
    title(['Velforce separability - ' num2str(velforce_sep)])
    figure(6)
    [velforce_sep,velforce_mdl] = test_sep(td,struct('signals',{{'glm_forcekin_pca'}},'mdl',actual_mdl,'do_plot',true));
    title(['Velforce w/ actual disc. separability - ' num2str(velforce_sep)])
    
    figure(7)
    [hand_elbow_sep,hand_elbow_mdl] = test_sep(td,struct('signals',{{'glm_hand_elbow_pca',1:4}},'do_plot',true));
    title(['hand_elbow separability - ' num2str(hand_elbow_sep)])
    figure(8)
    [hand_elbow_sep,hand_elbow_mdl] = test_sep(td,struct('signals',{{'glm_hand_elbow_pca'}},'mdl',actual_mdl,'do_plot',true));
    title(['hand_elbow w/ actual disc. separability - ' num2str(hand_elbow_sep)])

    % % then for directional separability/other view
    % [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    % [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    % signal_act = cat(1,td_act.glm_ext_pca);
    % signal_pas = cat(1,td_pas.glm_ext_pca);
    % % plot active as filled, passive as open
    % bump_colors = linspecer(4);
    % act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
    % pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;
    % subplot(2,3,2)
    % hold all
    % scatter3(signal_act(:,1),signal_act(:,2),signal_act(:,3),50,bump_colors(act_dir_idx,:),'filled')
    % scatter3(signal_pas(:,1),signal_pas(:,2),signal_pas(:,3),100,bump_colors(pas_dir_idx,:),'o','linewidth',2)
    % ylim = get(gca,'ylim');
    % zlim = get(gca,'zlim');
    % set(gca,'box','off','tickdir','out')
    % axis equal
    % axis off

%% get boostrapped separability values
    % get correlated noise models
    
    % bootstrap!
    n_boot = 1000;
    % use actual model for this
    bootsep_true = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'S1_FR_pca'}},'mdl',actual_mdl)),td');
    
    bootsep_ext = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'glm_ext_pca',1:4}})),td');
    % bootsep_ext_actual = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'glm_ext_pca'}},'mdl',actual_mdl)),td');
    
    bootsep_velforce = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'glm_forcekin_pca',1:4}})),td');
    % bootsep_velforce_actual = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'glm_forcekin_pca'}},'mdl',actual_mdl)),td');

    bootsep_handelbow = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'glm_hand_elbow_pca',1:4}})),td');
    % bootsep_handelbow_actual = bootstrp(n_boot,@(x) test_sep(x',struct('signals',{{'glm_hand_elbow_pca'}},'mdl',actual_mdl)),td');

    % get separability CIs
    sep_true = prctile(bootsep_true,[2.5 97.5]);
    
    sep_ext = prctile(bootsep_ext,[2.5 97.5]);
    % sep_ext_actual = prctile(bootsep_ext_actual,[2.5 97.5]);
    
    sep_velforce = prctile(bootsep_velforce,[2.5 97.5]);
    % sep_velforce_actual = prctile(bootsep_velforce_actual,[2.5 97.5]);

    sep_handelbow = prctile(bootsep_handelbow,[2.5 97.5]);
    % sep_handelbow_actual = prctile(bootsep_handelbow_actual,[2.5 97.5]);

    % plot separabilities with CIs
    figure('DefaultAxesFontSize',18)
    barh([mean(bootsep_true) mean(bootsep_ext)...
        mean(bootsep_velforce)...
        mean(bootsep_handelbow)],'edgecolor','none')
    % barh(fliplr([mean(bootsep_true) mean(bootsep_ext) mean(bootsep_ext_actual) ...
    %     mean(bootsep_velforce) mean(bootsep_velforce_actual) ...
    %     mean(bootsep_handelbow) mean(bootsep_handelbow_actual)]))
    hold on
    plot([sep_true' sep_ext'...
        sep_velforce'...
        sep_handelbow'],[1:4;1:4],'-k','linewidth',2)
    % plot(fliplr([sep_true' sep_ext' sep_ext_actual' ...
    %     sep_velforce' sep_velforce_actual' ...
    %     sep_handelbow' sep_handelbow_actual']),[1:7;1:7],'-k','linewidth',2)
    plot([0.5 0.5],[0 5],'--k','linewidth',2)
    plot([1 1],[0 5],'--k','linewidth',2)
    set(gca,'box','off','tickdir','out','yticklabel',{'Actual','Hand Kinematics','Handle Kinematics+Force',...
        'Hand/Elbow Kinematics'},...
        'xtick',[0 0.5 1],'xticklabel',{'','50%','100%'})
    axis ij

