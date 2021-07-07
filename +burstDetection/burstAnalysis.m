function Stats = burstAnalysis(DataOI_r,ala,features,max_bursts)
%%-------------------------------------------------
% performing analysis of burst features on trajectories in DataOI_r
%
%     DataOI_r (array shape n x m ) = dataset to plot
%        m categories
%        n trajectories (containing nan if not measured for this trajectory)
%
%     ala,features = return values from burstDetection(DataOI_r)
%
%     max_bursts (int) = [optional default=5] only first n bursts are taken into account
%
% returns 
%       burstStatistic object with fields:
%
%       Stats.maxBursts = maximal number of bursts (input argument); 
%
%     Statistics for bursts
%       Stats.height = height of bursts;
%       Stats.duration = duration of bursts;
%       Stats.interval = interval between bursts;
%       Stats.count = number of bursts;
%        
%
%     Statistics for population of trajectories:
%       Stats.popMean = mean of population for each time-point
%       Stats.popStd = std of population for each time-point
%       Stats.popStdVar = estimated variance of std (Jackknife resampled)
%       Stats.sqVar = quadratic variance of steps (measure of stochacisity)
%
% example:
%      [labels,bursts] = burstDetection(DataOI_r);
%      burstStats = burstAnalysis(DataOI_r,labels,bursts);
%
%%-------------------------------------------------
import burstDetection.burstStatistic
if ~exist('max_bursts','var')
    max_bursts = 5;
end
N = size(DataOI_r,2);

burst_volume = nan(N,max_bursts);
burst_height = nan(N,max_bursts);
burst_duration = nan(N,max_bursts);
burst_interval = nan(N,max_bursts);
burst_count = max(ala)';


for t = 1:N
    true;
    b = 1;
    for B = 1:size(features{t},1)
        hit = false;
        while b <= size(features{t},1) && ~hit
            if sum(ala(:,t)==b) > 0
                burst=DataOI_r(ala(:,t)==b,t);
                
                stop = features{t}(b,1);
                star = features{t}(b,2);
                dur = (stop - star)*5;
                hig = max(burst)-min(burst);
                
                if hig >= 0.03 && dur >= 15
                    hit = true;
                    if B <= max_bursts
                        burst_duration(t,B) = dur;
                        burst_height(t,B) = hig;
                        burst_volume(t,B) = sum(burst-min(burst))*5;


                        if B > 1
                            stop_last = features{t}(b-1,1);
                            burst_interval(t,B) = (star - stop_last)*5;
                        end
                    end
                    
                end

            end
            b = b + 1; % try next detected burst...
        end
        
        if ~hit
            burst_count(t) = burst_count(t) - 1;
        end
    end
end

N = size(DataOI_r,1);
n = round(N/max_bursts);
parts = [ 1:n:N ; [n:n:N,N]];
sq_var = nan(size(DataOI_r,2),max_bursts);
for i = 1:max_bursts
    seg = DataOI_r(parts(1,i):parts(2,i),:);
    sq_var = ( (seg(2:end,:)-seg(1:end-1,:)).^2 );
end

pop_mean = mean(DataOI_r,2);
pop_std = std(DataOI_r,[],2);

[T,N] = size(DataOI_r);
astar = zeros(T,N);
for t = 1:N
    astar(:,t) = std( DataOI_r(:,setdiff(1:N,t)),[],2 );
end
pop_std_var = sqrt( (N-1)/N*sum( (astar - mean(astar,2)).^2 ,2) );

%% wrap results in burstStatistic class object
Stats = burstStatistic();
 
Stats.maxBursts = max_bursts; 
Stats.height = burst_height;
Stats.duration = burst_duration;
Stats.interval = burst_interval;
Stats.count = burst_count;
Stats.popMean = pop_mean;
Stats.popStd = pop_std;
Stats.popStdVar = pop_std_var;
Stats.sqVar = sq_var;
 
end