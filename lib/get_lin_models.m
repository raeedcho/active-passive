function [lin_models] = get_lin_models(fr,thv,thf)
%% Find PDs of select neurons

% find each PD with circular mean
for neuron_ctr = 1:size(fr,2)
    dir_mat = [cos(thv) sin(thv) cos(thf) sin(thf)];
    fit_linear = LinearModel.fit(dir_mat,fr(:,neuron_ctr));

    lin_models{neuron_ctr} = fit_linear;
end