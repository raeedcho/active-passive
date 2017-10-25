%% set up TD
[~,td] = getTDidx(trial_data,'result','R');

td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));

[~,td_act] = getTDidx(td,'ctrHoldBump',false);
td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',15});
td_act = binTD(td_act,5);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);
td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_pas = binTD(td_pas,5);
td = cat(2,td_act,td_pas);

% scale spikes by binsize
for i = 1:length(td)
    td(i).S1_spikes = td(i).S1_spikes/td(i).bin_size;
end
[~,td_act] = getTDidx(td,'ctrHoldBump',false);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);

% opensim_idx = find(contains(td(1).opensim_names,'_vel'));
opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
signal_name = 'opensim';
num_boots = 1000;
opensimPDtable_act = getTDPDs(td_act,struct('out_signals',{{signal_name,opensim_idx}},'out_signal_names',{td(1).opensim_names(opensim_idx)},...
                                          'distribution','normal','move_corr','vel','num_boots',num_boots));
opensimPDtable_pas = getTDPDs(td_pas,struct('out_signals',{{signal_name,opensim_idx}},'out_signal_names',{td(1).opensim_names(opensim_idx)},...
                                          'distribution','normal','move_corr','vel','num_boots',num_boots));
spikesPDtable_act = getTDPDs(td_act,struct('out_signals',{{'S1_spikes'}},'out_signal_names',{td(1).S1_unit_guide},...
                                          'distribution','poisson','move_corr','vel','num_boots',num_boots));
spikesPDtable_pas = getTDPDs(td_pas,struct('out_signals',{{'S1_spikes'}},'out_signal_names',{td(1).S1_unit_guide},...
                                          'distribution','poisson','move_corr','vel','num_boots',num_boots));

% save PD tables
save('~/Projects/limblab/data-td/ForceKin/Han/20170203/pdtable50.mat','*PDtable*')

% plot tuning
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
subplot(1,2,2)
errorbar(opensimPDtable_act.velDir(isTuned),opensimPDtable_pas.velDir(isTuned),...
        minusPi2Pi(opensimPDtable_pas.velDir(isTuned)-opensimPDtable_pas.velDirCI(isTuned,1)),...
        minusPi2Pi(opensimPDtable_pas.velDirCI(isTuned,1)-opensimPDtable_pas.velDir(isTuned)),...
        minusPi2Pi(opensimPDtable_act.velDir(isTuned)-opensimPDtable_act.velDirCI(isTuned,1)),...
        minusPi2Pi(opensimPDtable_act.velDirCI(isTuned,1)-opensimPDtable_act.velDir(isTuned)),...
        'ro','linewidth',2)
hold on
plot([-pi pi],[-pi pi],'--k','linewidth',2)
xlabel('Active PD')
ylabel('Passive PD')
title('Muscle PDs')
set(get(gca,'ylabel'),'rotation',0,'horizontalalignment','right')
set(gca,'box','off','tickdir','out','xlim',[-pi pi],'ylim',[-pi pi],'xtick',[-pi pi],'ytick',[-pi pi],'xticklabel',{'-\pi','\pi'},'yticklabel',{'-\pi','\pi'})
axis equal

% plot neuron active tuning against passive tuning (same figure)
isTuned_params = struct('move_corr','vel','CIthresh',pi/3);
isTuned = checkIsTuned(spikesPDtable_act,isTuned_params)...
            & checkIsTuned(spikesPDtable_pas,isTuned_params);
subplot(1,2,1)
errorbar(spikesPDtable_act.velDir(isTuned),spikesPDtable_pas.velDir(isTuned),...
        (spikesPDtable_pas.velDir(isTuned)-spikesPDtable_pas.velDirCI(isTuned,1)),...
        (spikesPDtable_pas.velDirCI(isTuned,1)-spikesPDtable_pas.velDir(isTuned)),...
        (spikesPDtable_act.velDir(isTuned)-spikesPDtable_act.velDirCI(isTuned,1)),...
        (spikesPDtable_act.velDirCI(isTuned,1)-spikesPDtable_act.velDir(isTuned)),...
        'mo','linewidth',2)
hold on
plot([-pi pi],[-pi pi],'--k','linewidth',2)
xlabel('Active PD')
ylabel('Passive PD')
title('Neural PDs')
set(get(gca,'ylabel'),'rotation',0,'horizontalalignment','right')
set(gca,'box','off','tickdir','out','xlim',[-pi pi],'ylim',[-pi pi],'xtick',[-pi pi],'ytick',[-pi pi],'xticklabel',{'-\pi','\pi'},'yticklabel',{'-\pi','\pi'})
axis equal
