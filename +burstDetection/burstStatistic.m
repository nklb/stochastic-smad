classdef burstStatistic    
    % Container class for burst statistics
    properties
        %     Statistics for bursts
        %       height = height of bursts;
        %       duration = duration of bursts;
        %       interval = interval between bursts;
        %       count = number of bursts;
        %
        %
        %     Statistics for population of trajectories:
        %       popMean = mean of population for each time-point
        %       popStd = std of population for each time-point
        %       popStdVar = estimated variance of std (Jackknife resampled)
        %       sqVar = quadratic variance of steps (measure of 
        %               stochasticity)
        height = []
        duration = []
        interval = []
        count = []
        
        popMean = []
        popStd = []
        popStdVar = []
        sqVar = []
        
        maxBursts = []
    end
end

