import burstDetection.burstDetect;
import burstDetection.burstAnalysis;


data = load('data/Stimulation_100pM_2013.mat');
[label,bursts]=burstDetect(data.DataOI_r);
statistics = burstAnalysis(data.DataOI_r,label,bursts,4);

figure;
boxplot( statistics.height )
xlabel('# burst')
ylabel('Nuc/cyt SMAD2 ratio')
title('Height of burts')

import burstDetection.getRefStats
import burstDetection.objective_function
import forward.paths


ref = getRefStats(5,4,'2013');
simPaths = paths(data.DataOI_r,0:5:1440, [], 3, 5, []);
res = objective_function(simPaths,ref,4);