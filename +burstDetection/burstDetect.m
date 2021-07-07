function [labDetect,bursts] = burstDetect(trajectories,times,param)
%%-------------------------------------------------
% burst detection
%
% trajectories (array shape n x m) = set of m trajectories of size n
%
% times (array shape 1 x n) = timepoints of trajectory (modify if traj has
%                               6 min between each step)
%
% params (array of 8 floats) [optional default = optimized set] = set of 
%   parameters for burst detection
%   ! do not set params in usual application
%
% returns:
%    labDetect (array shape n x m) = labels: ala(t,k) > 0 if a burst was 
%       detected in trajectorie k at time t
%
%    bursts = set of n tuples describing positions (end, start) of burst k
%
%%-------------------------------------------------

if ~exist('times','var')
    times = 1:size(trajectories,1);
end
if ~exist('param','var')
    param = [49.8453   14.5377   10.1349   20.4543    3.8217   47.8792    0.0255    0.2776    0.0655   11.0694    0.7884];
end

filter_width_high = param(1);
ancor_weight_high = round(param(2));
filter_width_low = param(3);
ancor_weight_low =  round(param(4));
filter_width_final = param(5);
ancor_weight_final = round(param(6));
peak_min_prominence = param(7);
peak_min_dist = param(8);
peak_final_prominence = param(9);
peak_final_dist = param(10);
label_width_factor = param(11);


labDetect = zeros(size(trajectories,2),length(times));
bursts = cell(size(trajectories,2),1);

for k = 1:size(trajectories,2)
    
    smooth = gauss_est(trajectories(:,k),times,filter_width_low,ancor_weight_low);
    [~,posd,wd,hd] = findpeaks(smooth,times,'MinPeakProminence',peak_min_prominence,'MinPeakDistance',peak_min_dist);
    [~,total] = gauss_hills(times,posd,wd,hd);
    trend = gauss_est((smooth - total)',times,filter_width_high,ancor_weight_high);
    
    signal = trajectories(:,k)'-trend-min(trajectories(:,k)'-trend);
    smooth_signal = gauss_est(signal',times,filter_width_final,ancor_weight_final);
    [~,pos,w,~] = findpeaks(smooth_signal,times,'MinPeakProminence',peak_final_prominence,'MinPeakDistance',peak_final_dist);
        
    bursts{k} = zeros(length(pos),2);
    for p = 1:length(pos)
        support = times >= pos(p)-label_width_factor*w(p) & times <= pos(p)+label_width_factor*w(p);
        labDetect(k,support) = p;
        base = find(support);
        bursts{k}(p,:) = base([end,1]);%[pos(p),w(p),h(p)];
    end
    true;
end
labDetect = labDetect';
end

function [hills,total] = gauss_hills(times,pos,width,height)
% calculates a sum of gaussian hills to approximate a signal from parameters defined by findpeaks
    hills = zeros(length(pos),length(times));
    for i = 1:length(pos)
        hills(i,:) = height(i).*exp(-((times-pos(i))./width(i)).^2);
    end
    total = sum(hills,1);
    hills = hills';
end

function estimation = gauss_est(data_use,timepoints,ds,rpt)
%%-------------------------------------------------
% Gauss filter for input signal (used for burst detection )
% data_use (array): input signal
% timepoints (array): timepoints of input signal (defines distance)
% ds: filter size = width-parameter of gaussian kernel.
% rpt: initial value is repeated rpt times to emphasize the first datapoints.
%%-------------------------------------------------
    data_use = [data_use(1)*ones(rpt,1);data_use;data_use(end)*ones(rpt,1)];
    dt_i = timepoints(2)-timepoints(1);
    dt_e = timepoints(end)-timepoints(end-1);
    timepoints = [(timepoints(1)-rpt*dt_i):dt_i:(timepoints(1)-dt_i),timepoints,(timepoints(end)+dt_e):dt_e:(timepoints(end)+rpt*dt_e)];
    n_t = length(timepoints);
    estimation = zeros(1,n_t);
    for k = 1:n_t
        t_now = timepoints(k);
        
        contribution = pdf('Normal',timepoints,t_now,ds);
        weights = contribution'./sum(contribution);
        estimation(k) = sum(data_use.*weights);
    end
    estimation = estimation(rpt+1:end-rpt);
end
