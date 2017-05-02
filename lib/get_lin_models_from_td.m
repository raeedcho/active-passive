function lin_models = get_lin_models_from_td(td)
% might want to add a positional window if lin models aren't good enough
num_units = length(td(1).S1_unit_guide);
lin_models = cell(1,num_units);
for unitctr = 1:num_units
    thv = [];
    thf = [];
    fr = [];
    for trialctr = 1:length(td)
        % compose variables
        thv = [thv; atan2(td(trialctr).vel(:,2),td(trialctr).vel(:,1))];
        thf = [thf; atan2(td(trialctr).force(:,2),td(trialctr).force(:,1))];
        fr = [fr; td(trialctr).S1_spikes(:,unitctr)];
    end
    dir_mat = [cos(thv) sin(thv) cos(thf) sin(thf)];
    fr = fr/td(1).bin_size;
    fit_linear = LinearModel.fit(dir_mat,fr);

    lin_models{unitctr} = fit_linear;
end