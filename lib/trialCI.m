%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [avg_data, CIhigh_data, CIlow_data] = trialCI(trial_data, conditions, params)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [avg_data, CIhigh_data, CIlow_data] = trialCI(trial_data, conditions, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
do_stretch  =  false;
num_samp    =  1000;
nboot       =  1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
record_flag = true; % will add a flag field saying it's trial-averaged
if nargin > 2, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1, error('Conditions not provided as input.'); end
if ~iscell(conditions), conditions = {conditions}; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get list of time-varying signals that we will average over
time_vars = getTDfields(trial_data,'time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time warp each trial to the same number of points, if desired
if do_stretch
    for trial = 1:length(trial_data)
        for iVar = 1:length(time_vars)
            temp = trial_data(trial).(time_vars{iVar});
            trial_data(trial).(time_vars{iVar}) = interp1(1:size(temp,1),temp,linspace(1,size(temp,1),num_samp));
        end
    end
end

if length(unique(cellfun(@(x) size(x,1),{trial_data.pos}))) ~= 1
    error('Trials are not uniform length. Do this with time stretching option or trimTD.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get average data and condition indices
[avg_data,cond_idx] = trialAverage(trial_data,conditions,params);

% bootstrap averages
bootstat_td = bootstrp(nboot,@(x) trialAverage(x,conditions,params),trial_data);

% % go through time variables
% for i = 1:length(cond_idx)
%     % now loop along time signals to average
%     for v = 1:length(time_vars)
%         avg_data(i).(time_vars{v}) = mean(cat(3,trial_data(cond_idx{i}).(time_vars{v})),3);
%     end
%     if record_flag
%         avg_data(i).is_average = true;
%     end
%     % add idx fields if it wasn't time stretched
%     if ~do_stretch
%         fn_idx = getTDfields(trial_data,'idx');
%         for f = 1:length(fn_idx)
%             if length(unique([trial_data(cond_idx{i}).(fn_idx{f})])) == 1
%                 avg_data(i).(fn_idx{f}) = trial_data(cond_idx{i}(1)).(fn_idx{f});
%             end
%         end
%     end
% end
% 
% opensim_data = cat(3,temp_td.vel);
% opensimCI = prctile(opensim_data,[10 90],3);
% opensimCIlow = squeeze(opensimCI(:,:,1));
% opensimCIhigh = squeeze(opensimCI(:,:,2));
