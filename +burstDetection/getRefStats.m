function RefData = getRefStats(dose, numBursts, dataset)
% Provides reference statistics for the objective function 
%
% Pars:
% dose (integer between 0 and 6):   experiment indicator (0: 0 pM, 1: 1 pM,
%                                   2: 2.5 pM, 3: 5 pM, 4: 25 pM, 
%                                   5: 100 pM, 6: 2 x 5 pM)
% numBursts (integer):              number of bursts, for which height and
%                                   duration is analyzed
% dataset (string):                 part of the dataset:
%                                   '2013' = original dataset from 
%                                   Strasen paper
%                                   '2014' = equivalent experiments from 
%                                   2014
%                                   'combo' = 2013 + 2014 datasets combined
%
% Return values:
%   RefData (obj:burstStatistics):    reference statistics
% 
% example:
%   % Create a reference set for condition 5 (100 pM) using the first 4
%   % bursts from the 2013 original dataset only.
%   ref = getRefStats(5,4,'2013');
import burstDetection.getRootDir burstDetection.getRefPaths burstDetection.createRef

rootFolder = getRootDir();

if ~exist('dataset', 'var')
    dataset = '2014';
end

refStatFile = [rootFolder, 'data', filesep, 'refStat', filesep, ...
    'refStat', dataset, 'Dose', num2str(dose), 'Bursts', ...
    num2str(numBursts), '.mat'];

if exist(refStatFile, 'file')
    fprintf('Loading experimental reference data...\n');
    Loaded = load(refStatFile, 'RefData');
    RefData = Loaded.RefData;
else
    
    fprintf('Preparing experimental reference data...\n');
    refPaths = getRefPaths(dose, 1, dataset);
    RefData = createRef(refPaths, numBursts);
    save(refStatFile, 'RefData');
end