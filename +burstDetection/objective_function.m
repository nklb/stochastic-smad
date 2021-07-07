function res = objective_function(simPaths, RefStats, numPeaks, useControls)
%%-------------------------------------------------
% calculates measure of distance for burst properties
%
%   simPaths (paths object) =  object containing trajectories to compare
%                               to reference set
%
%   RefStats (burstStatistics object) = object containing cell statistics
%                                used as reference
%
%   numPeaks (int)  = number of bursts compared between simPaths and
%                       RefStats. Further bursts are ignored.
%
%   useControls (bool) = In simulated paths the data from control
%                       experiments is added to mimic measurement noise.
%                   
%                     true: Data from control experiments is added. 
%
%                     false: Data from control is not used for comparison.
%
% returns:
%    cost values  (array) = array of objective function components. 
%                            Dimension depends on number of bursts analyzed.
%
%       meaning of the array elements:
%
%           numPeaks dimensions: means of burst heights 
%           numPeaks dimensions: std of burst heights (weight = 0.5)
%
%           numPeaks dimensions: means of burst widths
%           numPeaks dimensions: td of burst widths (weight = 0.5)
%
%           1 dimension: mean of burst counts
%           1 dimension: std of burst counts (weight = 0.5*numPeaks)
%
%           1 dimension: L2 norm difference between population averages 
%               (weight = 10*numPeaks)
%           1 dimension: L2 norm difference between the standard-deviations 
%               of the population averages (weight = 10*numPeaks)
%
%           1 dimension: failed simulations (weight = 20*numPeaks);
% 
% example (assuming a function mySimuFunction that returns simulations as
%   path object):
% ref = getRefStats(5,4,'2013');
% simPaths = mySimuFunction();
%
% % compare simPaths to ref using 4 bursts only:
% res = objective_function(simPaths,ref,4)
%%-------------------------------------------------

import burstDetection.burstDetect burstDetection.burstStatistic burstDetection.burstAnalysis forward.paths

%% unwrap statistic of experimental data from burstStatistic class
height0 = RefStats.height;
width0 = RefStats.duration;
nburst0 = RefStats.count;
popmean0 = RefStats.popMean;
popstd0 = RefStats.popStd;

if ~exist('numPeaks','var')
    numPeaks = RefStats.maxBursts;
end

if ~exist('useControls','var')
    useControls = 1;
end

rel_fail = simPaths.getRelFiltered();
[~, genS] = simPaths.getPrepData(useControls);
[ala,feat] = burstDetect(genS);
SimStats = burstAnalysis(genS,ala,feat,numPeaks);

height = SimStats.height;
width = SimStats.duration;
nburst = SimStats.count;
popmean = SimStats.popMean;
popstd = SimStats.popStd;

[resHM, resHS] = featureDistMedStd(height0, height, numPeaks, 0);
[resDM, resDS] = featureDistMedStd(width0, width, numPeaks, 0);

[resNM, resNS] = featureDistMedStd(nburst0,nburst, 1, 1);

resM = featureDistL2(popmean0, popmean);

resS = featureDistL2(popstd0, popstd);

res = [resHM, 0.5 * resHS, resDM, 0.5 * resDS, numPeaks * resNM, numPeaks * 0.5 * resNS, numPeaks * 10 * resM, numPeaks * resS, numPeaks * 20 * rel_fail];%resI

res = res(:);

end

function [resMed, resStd] = featureDistMedStd(f0, f, numPeaks, use_mean)

if ~exist('use_mean','var')
    use_mean = false;
end
    
    resMed = zeros(1,numPeaks);
    resStd = zeros(1,numPeaks);
    
    for i = 1:numPeaks
        nanSim = isnan(f(:, i));
        nanDat = isnan(f0(:, i));
        
        rel_nan_sim = sum(nanSim)/length(nanSim);
        rel_nan_dat = sum(nanDat)/length(nanDat);
        rel_nan = subplus(rel_nan_sim - rel_nan_dat);
        
        if all(nanSim) || all(nanDat)
            resMed(i) = 1;
            resStd(i) = 1;
        else
            
            if use_mean
                medDat = mean(f0(~nanDat, i));
                medSim = mean(f(~nanSim, i));
            else
                medDat = quantile(f0(~nanDat, i), 0.5);
                medSim = quantile(f(~nanSim, i), 0.5);
            end
            
            resMed(i) = rel_nan + (1 - rel_nan) * abs(medDat - medSim)/(abs(medDat) + abs(medSim));
            
            stdDat = std(f0(~nanDat, i));
            stdSim = std(f(~nanSim, i));
            resStd(i) = rel_nan + (1 - rel_nan) * abs(stdDat - stdSim)/(abs(stdDat) + abs(stdSim));
        end
    end
end

function res = featureDistL2(f0, f)

res = norm(f0 - f, 1)/(norm(f0,1) + norm(f,1));

end