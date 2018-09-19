
%% Calculate tangling

% load some example data of neurons recorded from motor cortex and muscles
% data was recorded at 1kHz and includes a subset of the data from the two
% target pedaling task: few cycles of pedaling forward (condition 1) and a
% few cycles of pedaling backward (condition 2).
load sampleMtrCtx
load sampleEMG

[Q_m1, out_m1]   = tangleAnalysis(D_m1,  .001,'numPCs',8,'timeStep',5,'softenNorm', 5); % soft normalize neural data
[Q_emg, out_emg] = tangleAnalysis(D_emg, .001,'numPCs',8,'timeStep',5,'softenNorm', 0); % fully normalize EMG

% play around with the parameters (e.g. number of PCs) if you like. The
% only thing that cannot be changed is the sampling rate as it is intrinsic
% to the data.

%% Plot distributions of tangling (Q) across t1
histParams = {'BinWidth',100,'Normalization','cdf','DisplayStyle','stairs', 'BinLimits',[0 max(Q_emg)],'LineWidth',2};
figure; hold on;
h(1) = histogram(Q_m1, histParams{:}, 'EdgeColor','k');
h(2) = histogram(Q_emg, histParams{:}, 'EdgeColor','r');

axis([0 max(Q_emg) 0 1])
xlabel('tangling (Q)')
ylabel('probability')
legend(h,{'M1/PMd','EMG'},'Location','SouthEast')
set(gca,'FontSize',16)

%% Visualize tangling between pairs of time points (i.e. q(t1,t2) )
% 
tangle_visualize( out_m1 );% press any key to observe how tangling changes between pairs
% run this several times to see different t1 initializations
% Note that Q(t1) will be the maximum of q(t1,t2) across all t2

%% Same for EMG
tangle_visualize( out_emg ); 
