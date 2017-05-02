function plotDecomposedFR(fr,thv,thf,lin_models,savefigs)
% plots the firing rates at each useful reach of the OOR task and
% decomposes it into a velocity component, force component, and residual


% loop over units
for uid = 1:size(fr,2)
    

    % predict fr for partial models
    intercept = lin_models{uid}.Coefficients.Estimate(1);
    vel_coefs = lin_models{uid}.Coefficients.Estimate(2:3);
    force_coefs = lin_models{uid}.Coefficients.Estimate(4:5);
    
    vel_fr = [cos(thv) sin(thv)] * vel_coefs;
    force_fr = [cos(thf) sin(thf)] * force_coefs;
    
    resid_fr = fr(:,uid) - vel_fr - force_fr - intercept;
    
%     clim = [min(min(full_map)) max(max(full_map))];
    if savefigs
        h=figure('visible','off');
    else
        h = figure;
    end
    colormap(parula);
    subplot(221)
    tesselateFR(fr(:,uid),thv,thf)
    subplot(222)
    tesselateFR(vel_fr,thv,thf)
    subplot(223)
    tesselateFR(force_fr,thv,thf)
    subplot(224)
    tesselateFR(resid_fr,thv,thf)
    
    subplot(221)
    title(sprintf('Neuron %d', uid));
    
    if savefigs
        saveas(h,['C:\Users\rhc307\Projects\limblab\data-raeed\ForceKin\Han\20170203\Figures\Tesselate_fineLin\Neuron_' num2str(uid) '.png'])
    end
%     
%     close(h)
    
    disp(['Processed unit: ' num2str(uid)])
end % foreach unit

disp('Done!')

% figure
% imagesc(resid_map_total)
% colormap jet
