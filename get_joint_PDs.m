%% set up TD
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
td_act = binTD(td_act,15);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_pas = binTD(td_pas,15);
td = cat(2,td_act,td_pas);

% % get muscle velocity PCA
% PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel'))}},...
%                     'do_plot',true);
% [td,pca_info_vel] = getPCA(td,PCAparams_vel);
% % temporary hack to allow us to save into something useful
% for i=1:length(td)
%     td(i).opensim_muscVel_pca = td(i).opensim_pca;
% end

% opensim_idx = find(contains(td(1).opensim_names,'_vel'));
opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
signal_name = 'opensim';
num_boots = 10;
opensimPDtable_act = getTDPDs(td_act,struct('out_signals',{{signal_name,opensim_idx}},...
                                          'distribution','normal','move_corr','vel','num_boots',num_boots));
opensimPDtable_pas = getTDPDs(td_pas,struct('out_signals',{{signal_name,opensim_idx}},...
                                          'distribution','normal','move_corr','vel','num_boots',num_boots));
spikesPDtable_act = getTDPDs(td_act,struct('out_signals',{{'S1_spikes'}},...
                                          'distribution','normal','move_corr','vel','num_boots',num_boots));
spikesPDtable_pas = getTDPDs(td_pas,struct('out_signals',{{'S1_spikes'}},...
                                          'distribution','normal','move_corr','vel','num_boots',num_boots));

% plot tuning
maxRad = max(cat(1,opensimPDtable_act.velModdepth,opensimPDtable_pas.velModdepth));
colors = linspecer(length(opensim_idx));
h1 = figure;
polar(0,maxRad);
hold on
h2 = figure;
polar(0,maxRad);
hold on
for i = 1:length(opensim_idx)
    figure(h1)
    plotTuning([],opensimPDtable_act(i,:),[],opensimPDtable_act.velModdepth(i),colors(i,:),'-')
    figure(h2)
    plotTuning([],opensimPDtable_pas(i,:),[],opensimPDtable_pas.velModdepth(i),colors(i,:),'--')
end

% plot muscle active tuning against passive tuning
isTuned_params = struct('move_corr','vel','CIthresh',pi/4);
isTuned = checkIsTuned(opensimPDtable_act,isTuned_params)...
            & checkIsTuned(opensimPDtable_pas,isTuned_params);
figure
errorbar(opensimPDtable_act.velDir(isTuned),opensimPDtable_pas.velDir(isTuned),...
        minusPi2Pi(opensimPDtable_pas.velDir(isTuned)-opensimPDtable_pas.velDirCI(isTuned,1)),...
        minusPi2Pi(opensimPDtable_pas.velDirCI(isTuned,1)-opensimPDtable_pas.velDir(isTuned)),...
        minusPi2Pi(opensimPDtable_act.velDir(isTuned)-opensimPDtable_act.velDirCI(isTuned,1)),...
        minusPi2Pi(opensimPDtable_act.velDirCI(isTuned,1)-opensimPDtable_act.velDir(isTuned)),...
        'ro','linewidth',2)
hold on

% plot neuron active tuning against passive tuning (same figure)
isTuned_params = struct('move_corr','vel','CIthresh',pi/4);
isTuned = checkIsTuned(spikesPDtable_act,isTuned_params)...
            & checkIsTuned(spikesPDtable_pas,isTuned_params);
errorbar(spikesPDtable_act.velDir(isTuned),spikesPDtable_pas.velDir(isTuned),...
        (spikesPDtable_pas.velDir(isTuned)-spikesPDtable_pas.velDirCI(isTuned,1)),...
        (spikesPDtable_pas.velDirCI(isTuned,1)-spikesPDtable_pas.velDir(isTuned)),...
        (spikesPDtable_act.velDir(isTuned)-spikesPDtable_act.velDirCI(isTuned,1)),...
        (spikesPDtable_act.velDirCI(isTuned,1)-spikesPDtable_act.velDir(isTuned)),...
        'mo','linewidth',2)

plot([-pi pi],[-pi pi],'--k','linewidth',2)
xlabel('Active PD')
ylabel('Passive PD')
set(get(gca,'ylabel'),'rotation',0,'horizontalalignment','right')
set(gca,'box','off','tickdir','out','xlim',[-pi pi],'ylim',[-pi pi],'xtick',[-pi pi],'ytick',[-pi pi],'xticklabel',{'-\pi','\pi'},'yticklabel',{'-\pi','\pi'})
axis equal
