classdef paths
    % Class containing path data from simulation or experiments.
    properties (SetAccess = private)
        n                   % number of paths
        pathData            % raw path data (Nt x n array)
        t                   % time instances of data (1 x Nt array)
        
        fails = []          % indices of failed simulations (array)
        failMessages = []   % error messages for failed simulations (cell 
                            % array)
        
        complex             % indices of cell data including complex values
        negative            % indices of cell data including negative 
                            % values
        
        P                   % model parameters used in the simulations 
                            % (obj:modParameters)
        sigma               % noise parameters used in the simulations
                            % (102 x 2 array)
        dose                % experiment indicator (integer between 0 
                            % and 6)
        simDate = []        % time when the instance has been created
        
    end
    properties
        filterComplex = 1;  % if true, single cell data with complex values
                            % is filtered out
        filterNegative = 1; % if true, single cell data with negative
                            % values is filtered out
        useControls = 1;    % if true, zero stimulation control paths are
                            % added to the data
        info = [];          % additional information (string)
        experimental = 0;   % if true, instance includes experimental data
        
        TGFb = [];          % ligand data (Nt x 1 array). Not set by 
                            % default.
        dynVar = []         % data of dynamic concentrations (Nt x 23 
                            % array). Not set by default.
    end
    
    methods
        function PC = paths(pData, tIn, PIn, sigmaIn, doseIn, failsIn, ...
                failMessagesIn)
            % Constructor function
            %
            % Args:
            %   pData (N x m array):        nuc/cyt SMAD2 ratio of m cells 
            %                               at N time instances each
            %   tIn (1 x N array):          time instances 
            %                               (default: 0, 1, ..., N)
            %   Pin (obj:modParameters):    deterministic parameters used 
            %                               for the simulations
            %   sigmaIn (0, .., 6):         experimental case
            %   failsIn (m x 1 array):      indicator of failed simulations
            %                               (default: zero array)
            %   failMessagesIn 
            %   (1 x m cell array):         error messages for failed
            %                               path computations for debugging
            import forward.paths
            
            if ~exist('tIn', 'var')
                tIn =0:1:size(pData, 1)-1;
            end
            
            if numel(tIn) ~= size(pData, 1)
                error('Length of time vector must coincide with second dimension of path data.')
            end
            
            if exist('PIn', 'var')
                PC.P = PIn;
            end
            
            if exist('sigmaIn', 'var')
                PC.sigma = sigmaIn;
            end
            
            if exist('doseIn', 'var')
                PC.dose = doseIn;
            end
            
            PC.pathData = pData;
            PC.t = tIn;
            PC.n = size(pData, 2);
            
            
            if ~exist('failsIn', 'var')
                failsIn = zeros(PC.n, 1);
            end
            PC.fails = failsIn(:);
            PC.fails = find(PC.fails);
            
            if ~exist('failMessagesIn', 'var')
                failMessagesIn = cell(PC.n, 1);
            end
            failMessagesIn = failMessagesIn(:);
            PC.failMessages = failMessagesIn(~cellfun(@isempty, failMessagesIn));
            
            complexInd = ~all(isreal(pData), 1);
            PC.complex = find(complexInd).';
            
            negInd = ~all(pData >= 0, 1);
            PC.negative = find(negInd).';
            
            PC.simDate = datetime('now', 'TimeZone', 'local');
            
        end
        
        function out = allEmpty(PC)
            out = PC.n == numel(PC.getFilter());
        end
        
        function out = getRelFiltered(PC)
            out = numel(PC.getFilter())/PC.n;
        end
        
        function outPathData = getFiltered(PC)
            filter = getFilter(PC);
            
            valid = ~ismember(1:F.n, filter);
            
            outPathData = PC.pathData(:, valid);
            
        end
        
        function [tM, outPaths] = getPrepData(PC, controls)
            % Output path data after filtering and postprocessing
            %
            % Args:
            %   controls:                   if true zero stimulation 
            %                               control paths are added
            %                               (default: true for simulations,
            %                               false for experimental data)
            % Return values:
            % tM (1 x N array):             time instances of the path data
            % outPaths (m x N array):       path data  
            
            import burstDetection.getRootDir
            
            if ~exist('controls', 'var')
                controls = ~PC.experimental;
            end
            
            if PC.experimental
                if controls
                    error('Control paths cannot be used with experimental paths.');
                end
                tM = PC.t;
                outPaths = PC.pathData;
                return
            end
            
            filter = getFilter(PC);
            valid = ~ismember(1:PC.n, filter);
            rootFolder = getRootDir();
            
            if controls
                load([rootFolder, 'data', filesep, ...
                    'Centered_Controls.mat'], 'centered');
            end
            
            switch PC.dose
                case {1, 2, 3, 4, 5}
                    tM = 0:5:1440;
                case 6
                    tM = 0:5:1430;
                otherwise
                    error('Unknown time interval of reference data for this dose.');
            end
            
            if ~isempty(PC.pathData(:, valid))
                outPaths = interp1(PC.t, PC.pathData(:, valid), tM);
                if controls
                    idx = mod(0:PC.n-1, size(centered, 2)) + 1;
                    outPaths = outPaths + centered(1:size(outPaths, 1), ...
                        idx(valid));
                    if PC.filterNegative
                        outPaths = subplus(outPaths);
                    end
                end
            else
                outPaths = [];
            end
        end
        
        function av = popAv(PC, controls)
            % Output population average
            %
            % Args:
            %   controls:                   if true zero stimulation 
            %                               control paths are considered
            %                               (default: true for simulations,
            %                               false for experimental data)
            % Return values:
            % av (1 x N array):             population average of the
            %                               single cell simulations/data 
            
            if ~exist('controls', 'var')
                controls = ~PC.experimental;
            end
            
            [~, pData] = PC.getPrepData(controls);
            
            av = mean(pData, 2);
        end
        
        function F = visualize(paths, numShow, controls, popAv, expData, leg)
            % Visualizes the included paths
            %
            % Args:
            %   numShow (integer):  Number of paths to be shown. Chosen
            %                       randomly from all stored paths
            %                       (default: 20).
            %   controls (logical): adds control paths in the visualization
            %                       if true (default true for simulated
            %                       paths)
            %   popAv (logical):    adds the population average to the
            %                       figure if true (default: false)
            %   expData (logical):  adds the mean of the measurement data
            %                       for the corresponding experiment to the
            %                       figure if true (default: false)
            %   leg (logical):      adds a legend if true (default: false)
            %
            %   Returns:
            %   F (obj:Figure): Figure with paths.
            
            import burstDetection.getRefPaths
            if ~exist('numShow', 'var')
                numShow = 20;
            end
            
            if ~exist('controls', 'var')
                controls = ~paths.experimental;
            end
            
            if ~exist('popAv', 'var')
                popAv = false;
            end
            
            if ~exist('expData', 'var')
                expData = false;
            end
            
            if ~exist('leg', 'var')
                leg = false;
            end
            
            F = figure();
            
            [tAx, Paths] = paths.getPrepData(controls);
            numValPaths = size(Paths, 2);
            idx = randperm(numValPaths);
            Paths = Paths(:,idx(1:min(numShow, numValPaths)));
            
            tAx = tAx/60;
            if expData
                assert(~paths.experimental, 'This object already contains the experimental data.')
                Refpaths = getRefPaths(paths.dose, 0);
                RefMean = mean(Refpaths.pathData.');
                plot(tAx, RefMean,'k','LineWidth',2),hold on;
            end
            
            if popAv
                if paths.experimental
                    plot(tAx, paths.popAv(controls),'k','LineWidth',2),hold on;
                else
                    plot(tAx, paths.popAv(controls),'r','LineWidth',2),hold on;
                end
            end
            
            plot(tAx, Paths,':', 'Linewidth', 1.5);
            xlim([0 tAx(end)]);
            xlabel('time (h)');
            ylabel('nuc/cyt SMAD2')
            if leg
                if expData
                    if popAv
                        legend('average of data', 'average of simulation', 'simulation paths');
                    else
                        legend('population average');
                    end
                else
                    if popAv
                        legend('population average');
                    else
                        legend('simulation paths');
                    end
                end
            end
        end
        
    end
    methods(Access = private)
        function filter = getFilter(PC)
            % returns indices of paths to be filtered out
            filter = PC.fails;
            if PC.filterComplex
                filter = [filter; PC.complex];
            end
            if PC.filterNegative
                filter = [filter; PC.negative];
            end
            filter = unique(filter);
            
        end
    end
end