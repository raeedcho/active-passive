%% Use intended direction and force
[~,td] = getTDidx(trial_data,'result','R');
td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
td = trimTD(td,{'idx_movement_on',1},{'idx_movement_on',20});
% add firing rates rather than spike counts
td = addFiringRates(td,struct('array','S1'));
% td = binTD(td,20)

% fix 360 force direction
for i = 1:length(td)
    td(i).force_direction = mod(td(i).force_direction,360);
end

% get actual firing rates
for i = 1:length(td)
    td(i).S1_spikes = td(i).S1_spikes/td(i).bin_size;
end

% create models
[td,model_info] = getModel(td,struct('model_type','glm',...
    'model_name','S1_handle','in_signals',{{'vel'}},...
    'out_signals',{'S1_FR'}));
[td,model_info] = getModel(td,struct('model_type','glm',...
    'model_name','S1_velforce','in_signals',{{'vel';'force'}},...
    'out_signals',{'S1_FR'}));
opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
[td,model_info] = getModel(td,struct('model_type','glm',...
    'model_name','S1_muscle','in_signals',...
    {{'opensim',opensim_idx}},...
    'out_signals',{'S1_FR'}));

% trial average
td = trialAverage(td,{'target_direction','force_direction'});

% get target and load PDs
load_dirs = [td.force_direction]'/180*pi;
targ_dirs = [td.target_direction]'/180*pi;
for i = 1:size(td(1).S1_spikes,2)
    for condition = 1:length(td)
        spikes(condition) = mean(td(condition).S1_FR(:,i));
    end
    load_pds(i) = circ_mean(load_dirs,spikes')*180/pi;
    targ_pds(i) = circ_mean(targ_dirs,spikes')*180/pi;
end

for i = 1:size(td(1).S1_spikes,2)
    for condition = 1:length(td)
        force_idx = td(condition).force_direction/45 + 1;
        move_idx = td(condition).target_direction/45 + 1;
        
        tuning_mat(force_idx,move_idx,:) = td(condition).S1_spikes(:,i);
        tuning_mat_handle(force_idx,move_idx,:) = td(condition).glm_S1_handle(:,i);
        tuning_mat_velforce(force_idx,move_idx,:) = td(condition).glm_S1_velforce(:,i);
        tuning_mat_muscle(force_idx,move_idx,:) = td(condition).glm_S1_muscle(:,i);
    end
    tuning{i} = mean(tuning_mat,3);
    tuning_handle{i} = mean(tuning_mat_handle,3);
    tuning_velforce{i} = mean(tuning_mat_velforce,3);
    tuning_muscle{i} = mean(tuning_mat_muscle,3);
end

% interpolate tuning
tuning_recenter = zeros(size(tuning{1}));
[targ_mat,force_mat]=meshgrid(unique(load_dirs),unique(targ_dirs));
for i = 1:size(td(1).S1_spikes,2)
    tuning_recenter = tuning_recenter + interp2(1:8,1:8,tuning{i},minusPi2Pi(force_mat-load_pds(i))/(pi/4)+1,minusPi2Pi(targ_mat-targ_pds(i))/(pi/4)+1);
end
tuning_recenter = tuning_recenter/size(td(1).S1_spikes,2);

% h1 = figure;
% h2 = figure;
% h3 = figure;
% h4 = figure;
% n_rows = ceil(sqrt(size(td(1).S1_spikes,2)));
% for neuron_idx = 1:size(td(1).S1_spikes,2)
%     figure(h1)
%     subplot(n_rows,n_rows,neuron_idx)
%     imagesc(tuning{neuron_idx})
%     colorbar
%     clims = get(gca,'clim');
% 
%     figure(h2)
%     subplot(n_rows,n_rows,neuron_idx)
%     imagesc(tuning_handle{neuron_idx})
%     colorbar
%     caxis(clims)
% 
%     figure(h3)
%     subplot(n_rows,n_rows,neuron_idx)
%     imagesc(tuning_velforce{neuron_idx})
%     colorbar
%     caxis(clims)
% 
%     figure(h4)
%     subplot(n_rows,n_rows,neuron_idx)
%     imagesc(tuning_muscle{neuron_idx})
%     colorbar
%     caxis(clims)
% end

%% Try modeling with handle velocity

