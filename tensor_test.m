%% test out tensors on actpas data
    % First make tensor out of trial data
        [~,td] = getTDidx(trial_data,'result','R');
        td = removeBadTrials(td,struct('remove_nan_idx',false));
        td = getSpeed(td);
        td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
        td = smoothSignals(td,struct('signals','markers'));
        td = getDifferential(td,struct('signals','markers','alias','marker_vel'));
        % add firing rates rather than spike counts
        % td = softNormalize(td,'S1_spikes');
        td = smoothSignals(td,struct('signals',{{'S1_spikes'}},'calc_rate',true,'kernel_SD',0.05));
        td = sqrtTransform(td,'S1_spikes');

        % Remove unsorted channels
        keepers = (td(1).S1_unit_guide(:,2)~=0);
        % for Chips_20170913...
        if strcmpi(td(1).monkey,'Chips') && strcmpi(td(1).date,'09-13-2017')
            keepers(23) = false;
        end
        for trial = 1:length(td)
            td(trial).S1_unit_guide = td(trial).S1_unit_guide(keepers,:);
            td(trial).S1_spikes = td(trial).S1_spikes(:,keepers);
        end
        
        % get still handle data (no word for start of center hold)
        minCH = min(cat(1,td.ctrHold));
        bin_size = td(1).bin_size;
        still_bins = floor(minCH/bin_size);
        
        % Get td_act and td_pas
        % num_bins_before = floor(still_bins/2);
        num_bins_before = 15;
        num_bins_after = 30;
        
        % prep td_act
        [act_idx,td_act] = getTDidx(td,'ctrHoldBump',false);
        td_act = trimTD(td_act,{'idx_movement_on',-num_bins_before-5},{'idx_movement_on',num_bins_after-5});
        td(act_idx) = trimTD(td(act_idx),{'idx_movement_on',-num_bins_before},{'idx_movement_on',num_bins_after});
        
        % prep td_pas
        [pas_idx,td_pas] = getTDidx(td,'ctrHoldBump',true);
        td_pas = trimTD(td_pas,{'idx_bumpTime',-num_bins_before},{'idx_bumpTime',num_bins_after});
        td(pas_idx) = trimTD(td(pas_idx),{'idx_bumpTime',-num_bins_before},{'idx_bumpTime',num_bins_after});
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
        for trial = pas_idx
            assert(~isnan(td(trial).bumpDir),'Something went wrong...')
            td(trial).target_direction = td(trial).bumpDir*multiplier;
        end
        
        % concatenate tds
        td_cat = cat(2,td_act,td_pas);
        td_nocat = td;

        td = td_cat;

        % td = binTD(td,5);

    % make tensor
        num_factors = 10;

        % emg_data = cat(3,td.emg); % dimensions are timepoints x emg x trials
        % emg_data = emg_data-repmat(mean(mean(emg_data,3),1),size(emg_data,1),1,size(emg_data,3));
        % emg_tensor = tensor(emg_data);
        % M_emg = cp_als(emg_tensor,num_factors,'maxiters',100,'printitn',10);
        
        behave_data = cat(2,cat(3,td.marker_vel),cat(3,td.markers)); % dimensions are timepoints x markers x trials
        behave_data = behave_data-repmat(mean(mean(behave_data,3),1),size(behave_data,1),1,size(behave_data,3));
        behave_tensor = tensor(behave_data);
        M_behave = cp_als(behave_tensor,num_factors,'maxiters',200,'printitn',10);

        neural_data = cat(3,td.S1_spikes); % dimensions are timepoints x neurons x trials
        % neural_data = neural_data-repmat(mean(mean(neural_data,3),1),size(neural_data,1),1,size(neural_data,3));
        neural_tensor = tensor(neural_data);
        % M_neural = cp_als(neural_tensor,num_factors,'maxiters',300,'printitn',10);
        M_neural = cp_apr(neural_tensor,num_factors,'maxiters',300,'printitn',10);

    % color by direction
        num_colors = 4;
        dir_colors = linspecer(num_colors);
        trial_dirs = cat(1,td.target_direction);
        dir_idx = mod(round(trial_dirs/(2*pi/num_colors)),num_colors)+1;
        trial_colors = dir_colors(dir_idx,:);

    % plot tensor decomposition
        plotTensorDecomp(M_neural,struct('trial_colors',trial_colors,'bin_size',td(1).bin_size,'temporal_zero',num_bins_before/5+1))
