function out = getRefPaths(dose, rawData, dataset)
% Output reference data
% 
% Pars:
%   dose (int): selects the reference dataset 
%          1 = 1 pM
%          2 = 2.5 pM
%          3 = 5 pM
%          4 = 25 pM
%          5 = 100 pM
%          6 = 2x 5 pM (restimulation)
%
%  rawData (bool): return raw data array or path object
%          true: returns DataOI_r (raw data) array
%          false: returns paths object
%
%  dataset (string): which part of the dataset
%          2013 = original dataset from Strasen paper
%          2014 = equivalent experiments from 2014
%          combo = 2013 + 2014
%
% Return values:
%   out (obj:paths or array): reference data. Type depends on rawData
%                           parameter.
%%-------------------------------------------------

import burstDetection.getRootDir forward.paths

if ~exist('rawData', 'var')
    rawData = 0;
end

if ~exist('dataset', 'var') || strcmp(dataset, '2014')
    dataset = 14;
    dataset_info = ' from dataset Jette Brelin 2014.';
elseif strcmp(dataset, '2013')
    dataset = 13;
    dataset_info = ' from dataset Jette Brelin 2013.';
elseif strcmp(dataset, 'combo')
    dataset = 1314;
    dataset_info = ' from dataset Jette Brelin 2013 and 2014.';
else
    error('Dataset is not available!');
end

pathData = [];
dose_id = {'1pM','2_5pM','5pM','25pM','100pM', '2x2_5pM'};
if dose == 6 % && dataset ~= 14
    dataset = 20;
    L = load([getRootDir(), 'data', filesep, 'Stimulation_',dose_id{dose},'_2020.mat'], 'DataOI_r');
    pathData = [pathData,L.DataOI_r];
    dataset_info = ' from dataset Stefan Darmstadt 2020.';
    warning('There is only the restimulation data from 2020 available, which is loaded now.')
end
if dataset == 13 || dataset == 1314
    L = load([getRootDir(), 'data', filesep, 'Stimulation_',dose_id{dose},'_2013.mat'], 'DataOI_r');
    pathData = [pathData,L.DataOI_r];
end
if dataset == 14 || dataset == 1314
    L = load([getRootDir(), 'data', filesep, 'Stimulation_',dose_id{dose},'_2014.mat'], 'DataOI_r');
    pathData = [pathData,L.DataOI_r];
end
info = ['Measurement data for ',strrep(dose_id{dose},'_','.'),' Experiment, ',dataset_info];
t = 0:5:(size(pathData,1)-1)*5;

if rawData
    out = pathData;
else
    out = paths(pathData, t);
    out.info = info;
    out.experimental = true;
end
    
