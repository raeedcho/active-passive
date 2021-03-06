function [resid_lin,resid_non] = extract_resid(fr,thv,thf)
%% Find PDs of select neurons

% find each PD with circular mean
for neuron_ctr = 1:size(fr,2)
    dir_mat = [cos(thv) sin(thv) cos(thf) sin(thf)];
    fit_linear = LinearModel.fit(dir_mat,fr(:,neuron_ctr));

    coefs = fit_linear.Coefficients.Estimate;
    velPD = atan2(coefs(3),coefs(2));
    forcePD = atan2(coefs(5),coefs(4));
%     velPD = angle(fr(:,neuron_ctr)'*exp(1i*thv));
%     forcePD = angle(fr(:,neuron_ctr)'*exp(1i*thf));
%     
%     figure
%     polar(thv,fr(:,neuron_ctr),'o')
%     hold on
%     polar([velPD velPD], [0 10], 'r-')
    
    % fit other model
    new_mat = [cos(thv-velPD) cos(thf-forcePD) cos(thv-velPD).*cos(thf-forcePD)];
    new_fit = LinearModel.fit(new_mat,fr(:,neuron_ctr));
    
    % get residuals
    resid_lin(:,neuron_ctr) = fit_linear.Residuals.Raw;
    resid_non(:,neuron_ctr) = new_fit.Residuals.Raw;
end