function [separability,sep_train] = test_sep(td,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
use_trials      =  1:length(td);
signals         =  getTDfields(td,'spikes');
do_plot         =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters you can change that aren't described in header
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process and prepare inputs
signals = check_signals(td(1),signals);
if iscell(use_trials) % likely to be meta info
    use_trials = getTDidx(td,use_trials{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals = check_signals(td,signals);

% ideally, this would work like trialAverage, where the function would take in a condition
% But...I'm not sure how to deal with more than two values for the condition
[~,td_act] = getTDidx(td,'ctrHoldBump',false);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);

% clean nans out...?
nanners = isnan(cat(1,td_act.target_direction));
td_act = td_act(~nanners);

% plot active as filled, passive as open
bump_colors = linspecer(4);
signal_act = get_vars(td_act,signals);
signal_pas = get_vars(td_pas,signals);
act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;

% Find total separability
signal = cat(1,signal_act,signal_pas);
actpas = [ones(length(signal_act),1);zeros(length(signal_pas),1)];
[train_idx,test_idx] = crossvalind('LeaveMOut',length(actpas),floor(length(actpas)/10));
mdl = fitcdiscr(signal(train_idx,:),actpas(train_idx));
class = predict(mdl,signal(test_idx,:));
separability = sum(class == actpas(test_idx))/sum(test_idx);

class_train = predict(mdl,signal(train_idx,:));
sep_train = sum(class_train == actpas(train_idx))/sum(train_idx);

w = mdl.Sigma\diff(mdl.Mu)'
signal_sep = signal*w;

% get basis vector orthogonal to w for plotting
null_sep = null(w');
signal_null_sep = signal*null_sep;
[~,signal_null_sep_scores] = pca(signal_null_sep);

if do_plot
    figure
    hold all
    scatter3(signal_sep(actpas==1),signal_null_sep_scores(actpas==1,1),signal_null_sep_scores(actpas==1,2),50,bump_colors(act_dir_idx,:),'filled')
    scatter3(signal_sep(actpas==0),signal_null_sep_scores(actpas==0,1),signal_null_sep_scores(actpas==0,2),100,bump_colors(pas_dir_idx,:),'x','linewidth',3)
    % ylim = get(gca,'ylim');
    % plot([0 0],ylim,'--k','linewidth',2)
    set(gca,'box','off','tickdir','out')
    axis off
end
