function [compPaths] = compPathsVectorized(nPaths, P, sigma, dose, remNeg, ...
    sameSeed, dt, T, dispTime, newtonTol, newtonMaxIter, storeAllDyn)
% Time integration of the SMAD model using an implicit Taylor scheme
%
% Pars:
%   nPaths:                 number of paths to be computed
%   P (obj:modParameters):  model parameters as object of the class 
%                           modParameters
%   sigma (102x2 array):    noise parameters. First column: sigma, second
%                           column: theta, line index determines parameter.
%   dose:                   dose of TGFb (0: 0 pM, 1: 1 pM, 2: 2.5 pM,  
%                           3: 5 pM, 4: 25 pM, 5: 100 pM, 6: 2 x 5 pM)
%                           (default: 5)
%   remNeg:                 if true, paths that become negative will be 
%                           removed from the output
%   sameSeed:               if true, the seed for the random numbers is
%                           set to a fixed state before the simulation
%                           starts
%   dt:                     time increment for the scheme (default: 1)
%   T:                      final time of simutlation in minutes 
%                           (default: 1440)
%   disp_time:              if true, information about the computations
%                           is shown during simulation (default: false)
%   newtonTol:              tolerance of the newton method (default : 1e-6)
%   newtonMaxIter:          maximum newton iterations (default: 30)
%   storeAllDyn:            if true, the trajectories of all dynamic
%                           variables is stored in the output.
%                           Otherwise, only the nuc/cyt SMAD2 ratio is
%                           stored.
%
% Return values:
% compPaths (obj:paths):    Simulated paths as object of the class paths

import forward.expConditions forward.getInitialData ...
    forward.getSMADRatio forward.paths forward.SMADODE ...
    forward.sysJac forward.fastCIR

if ~exist('dose', 'var')
    dose = 5;
end

if ~exist('dt', 'var')
    dt = 1;
end

if ~exist('T', 'var')
    T = 1440;
end

if ~exist('dispTime', 'var')
    dispTime = 0;
end

if ~exist('newtonTol', 'var')
    newtonTol = 1e-6;
end

if ~exist('newtonMaxIter', 'var')
    newtonMaxIter = ceil(30*max(dt,1));
end

if ~exist('storeAllDyn', 'var')
    storeAllDyn = false;
end

if isempty(P.preDistribute)
    P.preDistribute = 0; % no predistribution as default
end



warnAlmostSing = warning('error', 'MATLAB:nearlySingularMatrix');
warnSing = warning('error', 'MATLAB:singularMatrix');

% load experimental conditions
E = expConditions(P, dose);

u0 = getInitialData(P);

%% actual simulation
if dispTime
    tic;
end

sInd = find(sigma(:, 1) ~= 0);
s0 = u0(sInd);
nStoch = length(sInd);

if sameSeed
    rng(87574, 'twister');
end

if P.preDistribute
    [~, sI, ~, ~] = fastCIR(s0, s0, sigma(sInd,:), nPaths, 1, 120);
    sStart = sI(:, end);
else
    sStart = repmat(s0, nPaths, 1);
end

%% compute stochastic parameters

if dose == 6
    t_part = [0, 2.5, 479, 483, T];
    dt_part = [dt/25, dt, dt/25, dt];
else
    t_part = [0, 2.5, T];
    dt_part = [dt/25, dt];
end

[t, s, dW, dZ] = fastCIR(sStart, s0, sigma(sInd,:), nPaths, dt_part(1), ...
    ceil(min(t_part(2), T)/dt_part(1)));
for i=2:length(dt_part)
    if T > t_part(i)
        [t2, s2, dW2, dZ2] = fastCIR(s(:,end), s0, sigma(sInd,:), nPaths, ...
            dt_part(i), ceil((min(t_part(i+1),T)-t_part(i))/dt_part(i)));
        t = [t(1:end-1); t2 + t_part(i)];
        s = [s(:, 1:end-1), s2];
        dW = [dW, dW2];
        dZ = [dZ, dZ2];
    end
end

Nt = length(t);
nDet = 23;
uD = zeros(nDet * nPaths, Nt);
uD(:, 1) = repmat(u0(1:nDet), nPaths, 1);

u = repmat(u0, nPaths, 1);

sig = repmat(sigma(sInd,1), nPaths, 1);
LPath = zeros(1, Nt);
free_lig = u0(5);
LPath(1) = free_lig;

%% construct help indices
K = 0:102:102*(nPaths-1);
dynInd = repmat(1:nDet, nPaths, 1) + diag(K) * ones(nPaths, nDet);
dynInd = sort(dynInd(:));

stochInd = repmat(sInd(:).', nPaths, 1) + diag(K) * ones(nPaths, length(sInd));
stochInd = sort(stochInd(:));

u(stochInd) = s(:, 1);

activeComp = false(102 * nPaths, 1);
activeDyn = false(nDet * nPaths, 1);

unFailedComp = true(102 * nPaths, 1);
unFailedDyn = true(nDet * nPaths, 1);
unFailedStoch = true(nStoch * nPaths, 1);

unFailedPaths = nPaths;
fail = false(nPaths, 1);

itInfo = zeros(1, newtonMaxIter);
fstatic = zeros(nDet * nPaths, 1);

for n = 1:Nt-1
    dt = t(n+1) - t(n);
    SM = SMADODE(t(n), u(unFailedComp), E, free_lig);
    if mod(n,3) == 1 || dt ~= t(n)-t(max(1, n-1))
        JcF = sysJac(t(n), u(unFailedComp), true);
    end
    fstatic(unFailedDyn) = uD(unFailedDyn, n) + 0.5 * dt * SM ...
        + JcF(:, stochInd(1:end-sum(fail)*nStoch)) ...
        * (sig(unFailedStoch) .* (dZ(unFailedStoch, n) ...
        - .5 * dt * dW(unFailedStoch, n)));
    u(stochInd) = s(:, n + 1);
    activeComp(:) = unFailedComp;
    activeDyn(:) = unFailedDyn;
    activePaths = unFailedPaths;
    for i = 1:newtonMaxIter
        if mod(i,3) == 1
            Jc = sysJac(t(n+1), u(activeComp), false);
        end
        JN = speye(nDet * activePaths) - 0.5 * dt * Jc;
        SM = SMADODE(t(n+1), u(activeComp), E, free_lig);
        rhs = u(dynInd(activeDyn))-fstatic(activeDyn) - 0.5 * dt * SM;
        
        if activePaths < 20
            try
                du = JN\rhs;
            catch
                itInfo(i) = activePaths;
                break
            end
            
        else
            prec = true;
            try
                [L,U] = ilu(JN);
            catch
                prec = false;
            end
            try
                if prec
                    [du,flag] = bicgstab(JN, rhs, 1e-12, 100, L, U);
                else
                    [du,flag] = bicgstab(JN, rhs, 1e-12, 100, L, U);
                end
            catch
                itInfo(i) = activePaths;
                break
            end
            if flag
                if dispTime
                    warning(['bicgstab failed, using LU decomposition', ...
                        'to solve the linear system.']);
                end
                try
                    du = JN\rhs;
                catch
                    itInfo(i) = activePaths;
                    break
                end
            end
        end
        
        u(dynInd(activeDyn)) = u(dynInd(activeDyn)) - du;
        unFailedCompInd = find(unFailedComp);
        free_lig = mean(u(unFailedCompInd(5:102:102*sum(unFailedPaths))));
        
        finished = sqrt(sum(reshape(rhs, nDet, activePaths).^2,1)) ...
            ./sqrt(sum(reshape(u(dynInd(activeDyn)), nDet, ...
            activePaths).^2,1)) < newtonTol;
        activePaths = activePaths - sum(finished);
        activeComp(activeComp) = logical(kron(~finished(:), ...
            true(102,1)));
        
        oaDyn = activeDyn;
        activeDyn(activeDyn) = logical(kron(~finished(:), ...
            ones(nDet,1)));
        itInfo(i) = activePaths;
        
        if activePaths && mod(i+1,3) ~= 1
            delComp = oaDyn & ~activeDyn;
            Jc(delComp(oaDyn), :) = [];
            Jc(:, delComp(oaDyn)) = [];
        elseif ~activePaths
            break
        end
    end
    uD(:, n+1) = u(dynInd);
    LPath(n+1) = free_lig;
    fail(activeDyn(1:nDet:end)) = true;
    
    ounFailedDyn = unFailedDyn;
    unFailedDyn(activeDyn) = false;
    delRow = ounFailedDyn & ~unFailedDyn;
    
    ounFailedComp = unFailedComp;
    unFailedComp(activeComp) = false;
    delCol = ounFailedComp & ~unFailedComp;
    
    JcF(delRow(ounFailedDyn), :) = [];
    JcF(:, delCol(ounFailedComp)) = [];
    
    unFailedStoch(logical(kron(fail(:), true(nStoch,1)))) = false;
    unFailedPaths = unFailedPaths - activePaths;
    
    if dispTime
        fprintf('t=%f, iterations: %s\n', t(n), num2str(itInfo(1:i)));
        itInfo(:) = 0;
    end
end

pathData = cell2mat(cellfun(@(Y)getSMADRatio(Y), ...
    mat2cell(uD, nDet * ones(nPaths, 1), size(uD,2)), ...
    'UniformOutput',false));

compPaths = paths(pathData.', t, P, sigma, dose, fail);
compPaths.filterNegative = remNeg;
compPaths.TGFb = LPath;


if dispTime
    toc;
end

if storeAllDyn
    compPaths.dynVar = [uD.', s.'];
end

warning(warnAlmostSing);
warning(warnSing);
