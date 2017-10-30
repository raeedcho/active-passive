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
opensimPDs{1} = getTDPDs(td_act,struct('out_signals',{{signal_name,opensim_idx}},'out_signal_names',{td(1).opensim_names(opensim_idx)},...
                                          'distribution','normal','move_corr','vel','num_boots',num_boots));
opensimPDs{2} = getTDPDs(td_pas,struct('out_signals',{{signal_name,opensim_idx}},'out_signal_names',{td(1).opensim_names(opensim_idx)},...
                                          'distribution','normal','move_corr','vel','num_boots',num_boots));
spikesPDs{1} = getTDPDs(td_act,struct('out_signals',{{'S1_spikes'}},'out_signal_names',{td(1).S1_unit_guide},...
                                          'distribution','poisson','move_corr','vel','num_boots',num_boots));
spikesPDs{2} = getTDPDs(td_pas,struct('out_signals',{{'S1_spikes'}},'out_signal_names',{td(1).S1_unit_guide},...
                                          'distribution','poisson','move_corr','vel','num_boots',num_boots));

% get tuning curves
[opensimCurves{1},bins] = getTuningCurves(td_act,struct('out_signals',{{signal_name,opensim_idx}},'out_signal_names',{td(1).opensim_names(opensim_idx)},'num_bins',4));
opensimCurves{2} = getTuningCurves(td_pas,struct('out_signals',{{signal_name,opensim_idx}},'out_signal_names',{td(1).opensim_names(opensim_idx)},'num_bins',4));
spikesCurves{1} = getTuningCurves(td_act,struct('out_signals',{{'S1_spikes'}},'out_signal_names',{td(1).S1_unit_guide},'num_bins',4));
spikesCurves{2} = getTuningCurves(td_pas,struct('out_signals',{{'S1_spikes'}},'out_signal_names',{td(1).S1_unit_guide},'num_bins',4));
% save PD tables
% save('~/Projects/limblab/data-td/ForceKin/Han/20170203/pdtable50.mat','*PDtable*')

% % plot tuning
% colors = linspecer(length(opensim_idx));
% h1 = figure;
% polar(0,maxRad);
% hold on
% h2 = figure;
% polar(0,maxRad);
% hold on
% for i = 1:length(opensim_idx)
%     figure(h1)
%     plotTuning([],opensimPDs{1}(i,:),[],opensimPDs{1}.velModdepth(i),colors(i,:),'-')
%     figure(h2)
%     plotTuning([],opensimPDs{2}(i,:),[],opensimPDs{2}.velModdepth(i),colors(i,:),'--')
% end

% plot muscle active tuning against passive tuning
isTuned_params = struct('move_corr','vel','CIthresh',pi/4);
isTuned = checkIsTuned(opensimPDs{1},isTuned_params)...
            & checkIsTuned(opensimPDs{2},isTuned_params);
figure
subplot(1,2,2)
errorbar(opensimPDs{1}.velDir(isTuned),opensimPDs{2}.velDir(isTuned),...
        minusPi2Pi(opensimPDs{2}.velDir(isTuned)-opensimPDs{2}.velDirCI(isTuned,1)),...
        minusPi2Pi(opensimPDs{2}.velDirCI(isTuned,1)-opensimPDs{2}.velDir(isTuned)),...
        minusPi2Pi(opensimPDs{1}.velDir(isTuned)-opensimPDs{1}.velDirCI(isTuned,1)),...
        minusPi2Pi(opensimPDs{1}.velDirCI(isTuned,1)-opensimPDs{1}.velDir(isTuned)),...
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
isTuned = checkIsTuned(spikesPDs{1},isTuned_params)...
            & checkIsTuned(spikesPDs{2},isTuned_params);
subplot(1,2,1)
errorbar(spikesPDs{1}.velDir(isTuned),spikesPDs{2}.velDir(isTuned),...
        (spikesPDs{2}.velDir(isTuned)-spikesPDs{2}.velDirCI(isTuned,1)),...
        (spikesPDs{2}.velDirCI(isTuned,1)-spikesPDs{2}.velDir(isTuned)),...
        (spikesPDs{1}.velDir(isTuned)-spikesPDs{1}.velDirCI(isTuned,1)),...
        (spikesPDs{1}.velDirCI(isTuned,1)-spikesPDs{1}.velDir(isTuned)),...
        'mo','linewidth',2)
hold on
plot([-pi pi],[-pi pi],'--k','linewidth',2)
xlabel('Active PD')
ylabel('Passive PD')
title('Neural PDs')
set(get(gca,'ylabel'),'rotation',0,'horizontalalignment','right')
set(gca,'box','off','tickdir','out','xlim',[-pi pi],'ylim',[-pi pi],'xtick',[-pi pi],'ytick',[-pi pi],'xticklabel',{'-\pi','\pi'},'yticklabel',{'-\pi','\pi'})
axis equal

% plot all tuning
figure
plotAllTuning(opensimCurves,opensimPDs,bins,find(isTuned))
figure
plotAllTuning(spikesCurves,spikesPDs,bins,find(isTuned))

