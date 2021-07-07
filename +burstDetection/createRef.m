function RefData = createRef(DataOI_r,max_bursts,file,enforce_recalc)
%%-------------------------------------------------
% applies burst detection and burst analysis to reference set (experimental dataset)
% and saves burst properties so that they can be used in objective_function
%
%    DataOI_r (array shape n x m) = set of m trajectories of size n
%
%    max_bursts (int) = number of bursts analyzed
%
%    file (string) = file name to load/save result of burstDetection
%    from/to
%
%    enforce_recalc (bool) = if true: ignore existing file, recompute
%       bursts using burstDetection and overwrite file
%
% returns:
%     RefData (obj:burstStatistic) = object with fields for different statistics
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
%%-------------------------------------------------
   import burstDetection.burstDetect burstDetection.burstAnalysis burstDetection.burstStatistic
    
    if ~exist('enforce_recalc','var')
        enforce_recalc = false;
    end
    if ~exist('file','var') || isempty( file ) || enforce_recalc
        [ala0,feat0]=burstDetect(DataOI_r);
        if exist('file','var') && ~isempty( file )
            save(file,'ala0','feat0');
        end
    else
        load(file,'ala0','feat0');
    end
    
    if ~exist('max_bursts','var')
        max_bursts = 5;
    end
    
    RefData = burstAnalysis(DataOI_r,ala0,feat0,max_bursts);
    
end