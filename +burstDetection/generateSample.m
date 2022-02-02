function [hills_added,label] = generateSample(sample_size,diff_crit,noise_scale,x_min,x_max,max_peaks,width_min,width_range,height_min,height_range,show_plot)
%%-------------------------------------------------
% generates trajectories with artificial bursts to use as trainingsset for
% burst detection and analyzis
%
%     sample_size (int) = number of trajectories generated
% 
%     for remaining parameters see short descriptions in comments where
%     their defaults get initialized.
%
% returns:
% dGN (array shape n x m) = set of m trajectories of length n.
%
% label (array shape n x m) = label(t,k) is > 0 if trajectorie k shows a
%                   burst at time point t. this is equivalent to ala from
%                   featP53.
%
% example:
% [samp,label] = generateSample(100);
% [ala,features] = featP53(samp);
% bursts_found = ala > 0 & label > 0;  % contains true positives of featP53
% bursts_missed = ala == 0 & label > 0; % contains false negatives of featP53 
% 
%%-------------------------------------------------
if ~exist('sample_size','var')
    sample_size = 2000; % size of sample
end
if ~exist('diff_crit','var')
    diff_crit = 4.5; % minimal distance of peaks normalized to peaks width
end
if ~exist('noise_scale','var')
    noise_scale = 0.015; % noise 2
end
if ~exist('max_peaks','var')
    max_peaks = 16;  % maximal number of bursts
end
if ~exist('x_min','var')
    x_min = 55; % minimum x for peak position to avoid burst in signal-peak
end
if ~exist('x_max','var')   
    x_max = 260; % maximum x for peak position 
end
if ~exist('width_min','var')
    width_min = 4; % minimal width of peaks
end
if ~exist('width_range','var')
    width_range = 6; % maximal_width = width_min + width_range 
end
if ~exist('height_min','var')
    height_min = 3;  % minimal height of peaks
end
if ~exist('height_range','var')
    height_range = 2; % maximal_height = height_min + height_range 
end
if ~exist('show_plot','var')
    show_plot = false; % bool: show a set of generated trajectories for visual quality check
end

load('data/Stimulation_100pM_2013.mat');
population_mean=mean(DataOI_r,2);
len = length(population_mean);
label = zeros(len,sample_size);
hills_added = zeros(len,sample_size);

[classes,archetypes]=kmeans(DataOI_r',3);
prop_archetype = zeros(1,size(archetypes,1));
for i = 1:size(archetypes,1)
    prop_archetype(i)  = sum(classes==i)/length(classes);
end

% figure,plot(arctypes');
% legend(num2str(sets(1)),num2str(sets(2)),num2str(sets(3)))


for k = 1:sample_size
    take_archetype = randsample(1:size(archetypes,1),1,true,prop_archetype);
    base_trajectory = archetypes(take_archetype,:)';
    noise=cumsum(randn(length(base_trajectory),1).*noise_scale);
    trajectory = base_trajectory + noise;
    
    features = [];
    n_peaks = randi([1,max_peaks]);
    for p = 1:n_peaks
        % generate properties of next burst
        pos = randi([x_min,x_max]);
        width = width_range*rand()+width_min;
        height = height_range*rand()+height_min;
    
        % to avoid peaks that are too close to each other to be discriminated,
        % a measure of discriminability from previously added peaks is evaluated
        diff_to_next_peak = [];
        for existing_burst = 1:size(features,1)
            pos_old = features(existing_burst,1);
            width_old = features(existing_burst,2);
            diff_to_next_peak(end+1) = abs((pos-pos_old)/sqrt(width*width_old));
        end
        
        % add next peak only if discriminatable and save features
        if isempty(diff_to_next_peak) || min(diff_to_next_peak) > diff_crit
            burst = (pdf('Normal',1:len,pos,width)').*height;
            trajectory = trajectory+burst;
            features(end+1,:) =  [pos,width];
            
            start = max([1,ceil( pos - 2*width )]);
            stop  = min([len,floor( pos + 2*width )]);
            label(start:stop,k) = size(features,1)+1;
        end
    end
    label(7:40,k) = 1; % label initial peak from base signal
    hills_added(:,k) = trajectory;
    
    %% plot examples
    if show_plot && k <= 9
        figure(3);
        if k == 1
            clf;
        end
        subplot(3,3,k);
        dx = 0:5:1390;
        plot(dx,trajectory,'b-');hold on
        % plot labeled hills
        for l = 1:(length(features)+1)
            take_this = label(:,k)==l;
            plot(dx(take_this),trajectory(take_this),'r-','LineWidth',2);
        end
        title(sprintf('%d bursts',max(label(:,k))));
        if ~isempty(intersect([1,4,7],k))
            ylabel('nuc/cyt SMAD2 ratio');
        end
        if ~isempty(intersect([7,8,9],k))
            xlabel('time / min');
        end
        
    end
end

end

