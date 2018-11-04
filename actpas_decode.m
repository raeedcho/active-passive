%% setup
    [~,td] = getTDidx(trial_data,'result','R');
    
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
    
    % bin trial data, then shift, then trim
    td = binTD(td,5);
    td = dupeAndShift(td,'S1_spikes',-3);

    % Get td_act and td_pas
    num_bins_before = 3;
    num_bins_after = 5;
    
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

    % clean up
    clearvars -except trial_data td

%% Crossval fit and predict decoders
    num_repeats = 20;
    num_folds = 5;

    pos_model_params = struct(...
        'model_type','linmodel',...
        'model_name','pos',...
        'in_signals',{{'S1_spikes','all';'S1_spikes_shift','all'}},...
        'out_signals',{{'pos','all'}});
    vel_model_params = struct(...
        'model_type','linmodel',...
        'model_name','vel',...
        'in_signals',{{'S1_spikes','all';'S1_spikes_shift','all'}},...
        'out_signals',{{'vel','all'}});

    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    minsize = min(length(td_act),length(td_pas));
    td_act = td_act(1:minsize);
    td_pas = td_pas(1:minsize);

    td_test = cell(num_repeats,num_folds);
    pos_eval = cell(num_repeats,num_folds);
    vel_eval = cell(num_repeats,num_folds);
    repeat_tic = tic;
    for repeatnum = 1:num_repeats
        xval_idx = crossvalind('kfold',minsize,num_folds);
        fold_tic = tic;
        for foldnum = 1:num_folds
            test_idx = (xval_idx==foldnum);
            td_test{repeatnum,foldnum} = cat(2,td_act(test_idx),td_pas(xval_idx));
            td_train = cat(2,td_act(~test_idx),td_pas(~test_idx));

            [~,pos_model] = getModel(td_train,pos_model_params);
            [~,vel_model] = getModel(td_train,vel_model_params);
            
            td_test{repeatnum,foldnum} = getModel(td_test{repeatnum,foldnum},pos_model);
            td_test{repeatnum,foldnum} = getModel(td_test{repeatnum,foldnum},vel_model);

            pos_model.eval_metric = 'vaf';
            vel_model.eval_metric = 'vaf';
            pos_model.num_boots = 1;
            vel_model.num_boots = 1;
            pos_eval{repeatnum,foldnum} = evalModel(td_test{repeatnum,foldnum},pos_model);
            vel_eval{repeatnum,foldnum} = evalModel(td_test{repeatnum,foldnum},vel_model);

            fprintf('\tEvaluated fold %d of %d at time %f\n',foldnum,num_folds,toc(fold_tic))
        end
        fprintf('Evaluated repeat %d of %d at time %f\n',repeatnum,num_repeats,toc(repeat_tic))
    end

    %% make plot
    td_plot = horzcat(td_test{end,:});
    true_pos = getSig(td_plot,{'pos',1:2});
    true_vel = getSig(td_plot,{'vel',1:2});
    pred_pos = getSig(td_plot,{'linmodel_pos',1:2});
    pred_vel = getSig(td_plot,{'linmodel_vel',1:2});

    figure('defaultaxesfontsize',18)
    ax(1) = subplot(2,1,1);
    scatter(true_pos(:,1),pred_pos(:,1),[],'k','filled')
    hold on
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    plot(xlim,[0 0],'-k','linewidth',2)
    plot([0 0],ylim,'-k','linewidth',2)
    plot(xlim,ylim,'--k','linewidth',2)
    title('Predicted vs Actual Position')
    ylabel('Predicted X-coordinate (cm)');
    xlabel('True X-coordinate (cm)');
    set(gca,'tickdir','out','box','off')
    ax(2) = subplot(2,1,2);
    scatter(true_pos(:,2),pred_pos(:,2),[],'k','filled')
    hold on
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    plot(xlim,[0 0],'-k','linewidth',2)
    plot([0 0],ylim,'-k','linewidth',2)
    plot(xlim,ylim,'--k','linewidth',2)
    title('Predicted vs Actual Position')
    ylabel('Predicted Y-coordinate (cm)');
    xlabel('True Y-coordinate (cm)');
    set(gca,'tickdir','out','box','off')

    figure('defaultaxesfontsize',18)
    subplot(2,1,1)
    plot([0 1],[0 0],'-k','linewidth',2)
    hold on
    plot([0 0],[0 1],'-k','linewidth',2)
    plot([0 1],[0 1],'--k','linewidth',2)
    set(gca,'tickdir','out','box','off')
    subplot(2,1,2)
    plot([0 1],[0 0],'-k','linewidth',2)
    hold on
    plot([0 0],[0 1],'-k','linewidth',2)
    plot([0 1],[0 1],'--k','linewidth',2)
    set(gca,'tickdir','out','box','off')
    for repeatnum = 1:num_repeats
        for foldnum = 1:num_folds
            subplot(2,1,1)
            scatter(pos_eval{repeatnum,foldnum}(1),vel_eval{repeatnum,foldnum}(1),[],'k','filled')
            subplot(2,1,2)
            scatter(pos_eval{repeatnum,foldnum}(2),vel_eval{repeatnum,foldnum}(2),[],'k','filled')
        end
    end

