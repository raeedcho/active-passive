%% Make rasters
    %% Set up trial data
        [~,td] = getTDidx(trial_data,'result','R');

        % Remove unsorted channels
        keepers = (td(1).S1_unit_guide(:,2)~=0);
        for trial = 1:length(td)
            td(trial).S1_unit_guide = td(trial).S1_unit_guide(keepers,:);
            td(trial).S1_spikes = td(trial).S1_spikes(:,keepers);
        end

        % remove low firing neurons
        td = removeBadNeurons(td,struct('min_fr',0.1));
        
        td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
        td = smoothSignals(td,struct('signals','markers'));
        td = getDifferential(td,struct('signals','markers','alias','marker_vel'));
        % add firing rates rather than spike counts
        td = addFiringRates(td,struct('array','S1'));
        
        % split into act and pas, then trim to length
        num_bins_before = 15;
        num_bins_after = 30;

        [~,td_act] = getTDidx(td,'ctrHoldBump',false);
        % clean nans out...?
        nanners = isnan(cat(1,td_act.target_direction));
        td_act = td_act(~nanners);
        td_act = trimTD(td_act,{'idx_movement_on',-num_bins_before},{'idx_movement_on',num_bins_after-1});

        [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
        td_pas = trimTD(td_pas,{'idx_bumpTime',-num_bins_before},{'idx_bumpTime',num_bins_after-1});
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

        % even out sizes
        minsize = min(length(td_act),length(td_pas));
        td_act = td_act(1:minsize);
        td_pas = td_pas(1:minsize);

    %% Fit and predict models...
        model_aliases = {'ext','handelbow'};
        model_params = cell(length(model_aliases),1);
        for modelnum = 1:length(model_aliases)
            switch model_aliases{modelnum}
            case 'ext'
                markername = 'Marker_1';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td(1).marker_names);
                assert(all(point_exists),'Hand marker does not exist?')
                model_params{modelnum} = struct(...
                    'model_type','glm',...
                    'model_name',model_aliases{modelnum},...
                    'in_signals',{{'markers',marker_hand_idx;'marker_vel',marker_hand_idx}},...
                    'out_signals','S1_FR');
            case 'handelbow'
                markername = 'Marker_1';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td(1).marker_names);
                assert(all(point_exists),'Hand marker does not exist?')

                markername = 'Pronation_Pt1';
                [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td(1).marker_names);
                assert(all(point_exists),'Elbow marker does not exist?')
                model_params{modelnum} = struct(...
                    'model_type','glm',...
                    'model_name',model_aliases{modelnum},...
                    'in_signals',{{'markers',[marker_hand_idx marker_elbow_idx];'marker_vel',[marker_hand_idx marker_elbow_idx]}},...
                    'out_signals','S1_FR');
            end
        end

        % cross-validate to get predictions
        num_folds = 5;
        xval_idx = crossvalind('kfold',minsize,num_folds);
        td_test = cell(1,num_folds);
        foldtic = tic;
        for foldnum = 1:num_folds
            test_idx = (xval_idx == foldnum);
            train_idx = ~test_idx;

            td_test{foldnum} = cat(2,td_act(test_idx),td_pas(test_idx));
            td_train = cat(2,td_act(train_idx),td_pas(train_idx));

            for modelnum = 1:length(model_aliases)
                [~,model_info] = getModel(td_train,model_params{modelnum});

                td_test{foldnum} = getModel(td_test{foldnum},model_info);
            end
            fprintf('Evaluated fold %d at time %f\n',foldnum,toc(foldtic))
        end
        td = horzcat(td_test{:});

        % normalize and split again
        % td = smoothSignals(td,struct('signals','S1_spikes','calc_rate',true));
        % td = sqrtTransform(td,'S1_spikes');
        % td = sqrtTransform(td,'S1_FR');
        td = softNormalize(td,struct('signals','S1_FR','alpha',1));
        for modelnum = 1:length(model_aliases)
            td = softNormalize(td,struct('signals',strcat('glm_',model_aliases{modelnum}),'alpha',1));
        end
        % td = zscoreSignals(td,struct('signals','S1_FR'));
        [~,td_act] = getTDidx(td,'ctrHoldBump',false);
        [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

    %% Average over act/pas and target direction
        td_act_avg = trialAverage(td_act,'target_direction');
        td_pas_avg = trialAverage(td_pas,'target_direction');

    %% Plot out average "raster" for direction of movement
        % pick which raster to plot
        % models_to_plot = {'S1_FR','glm_ext','glm_handelbow'};
        models_to_plot = {'S1_FR'};
        num_dirs = 4;

        cm_viridis = viridis(200);

        % make plots
        % h=cell(num_dirs,1);
        % for trialnum = 1:num_dirs
        %     h{trialnum} = figure('defaultaxesfontsize',18);
        % end
        figure('defaultaxesfontsize',18)
        binvec = (0:num_bins_before:num_bins_before*3);
        timevec = (binvec-num_bins_before)*td(1).bin_size*1000;
        for modelnum = 1:length(models_to_plot)
            fullraster = [[td_act_avg.(models_to_plot{modelnum})],...
                [td_pas_avg.(models_to_plot{modelnum})]];
            clim = [min(min(fullraster)) max(max(fullraster))];
            for trialnum = 1:num_dirs
                % figure(h{trialnum})
                % plot active
                raster = td_act_avg(trialnum).(models_to_plot{modelnum})';
                subplot(num_dirs,2,2*(trialnum-1)+1)
                imagesc(raster,clim);
                hold on
                plot([1 1]*num_bins_before,[0 size(raster,1)+1],'--w','linewidth',3)
                plot([1 1]*(num_bins_before+15),[0 size(raster,1)+1],'--w','linewidth',3)
                hold off
                set(gca,...
                    'box','off',...
                    'tickdir','out',...
                    'xtick',binvec,...
                    'xticklabel',timevec,...
                    'ytick',[]);
                % xlabel('Time after movement onset (ms)')
                % ylabel('Neuron')
                ylabel(sprintf('Target Direction: %f',td_act_avg(trialnum).target_direction))

                % plot passive
                raster = td_pas_avg(trialnum).(models_to_plot{modelnum})';
                subplot(num_dirs,2,2*(trialnum-1)+2)
                imagesc(raster,clim);
                hold on
                plot([1 1]*num_bins_before,[0 size(raster,1)+1],'--w','linewidth',3)
                plot([1 1]*(num_bins_before+15),[0 size(raster,1)+1],'--w','linewidth',3)
                hold off
                set(gca,...
                    'box','off',...
                    'tickdir','out',...
                    'xtick',binvec,...
                    'xticklabel',timevec,...
                    'ytick',[]);
                % xlabel('Time after bump onset (ms)')
                % ylabel('Neuron')
                ylabel(sprintf('Target Direction: %f',td_pas_avg(trialnum).target_direction))

                colormap(cm_viridis)
            end
        end

