function td = predict_td_from_models(td,lin_models)

num_units = length(lin_models);
% actpas_unit_idx = find(td(1).S1_unit_guide(:,2)~=0);
if num_units ~= length(td(1).S1_unit_guide)
    error('wrong number of units')
end

for trialctr = 1:length(td)
    % clear the current S1 spikes
    td(trialctr).S1_spikes = zeros(size(td(trialctr).S1_spikes));
    for unitctr = 1:num_units
        % predict firing rates
        thv = atan2(td(trialctr).vel(:,2),td(trialctr).vel(:,1));
        thf = atan2(td(trialctr).force(:,2),td(trialctr).force(:,1));
        dir_mat = [cos(thv) sin(thv) cos(thf) sin(thf)];
        pred_fr = lin_models{unitctr}.predict(dir_mat);
        
        % rescale for spike counts
        pred_spikes = pred_fr*td(trialctr).bin_size;
        pred_spikes(pred_spikes<0) = 0;
        
        % assign predictions to td
        td(trialctr).S1_spikes(:,unitctr) = pred_spikes;
    end
end