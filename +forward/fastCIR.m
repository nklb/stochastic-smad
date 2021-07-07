function [t, u, dW, dZ] = fastCIR(u0, uA, sigma, nPaths, dt, Nt)
% Computation of CIR processes. Used in the path computation
%
% Pars:
%   u0 (mx1 array):         initial value of the processes
%   uA (mx1 array):         attracting state of the processes (in most 
%                           cases equal to u0)            
%   sigma (mx2 array):      noise parameters. First column: sigma, second
%                           column: theta.
%   nPaths:                 number of processes to be computed per noise
%                           parameter set
%   dt:                     time increment in the process simulation
%   Nt:                     number of time steps to be simulated
%
% Return values:
% t (1x[Nt+1] array):           time instances of the processes
% u ([m*nPaths]x[Nt+1] array):  CIR process values
% dW ([m*nPaths]xNt array):     first set of random numbers used in the 
%                               simulation
% dW ([m*nPaths]xNt array):     second set of random numbers used in the 
%                               simulation

t = (0:Nt).' * dt;

nStoch = size(sigma,1);
u = zeros(nStoch * nPaths, Nt+1);

if length(u0) == nStoch
    u(:,1) = repmat(u0, nPaths, 1);
elseif length(u0) == nStoch * nPaths
    u(:,1) = u0;
else
    error('Invalid length of u0.')
end
uAtr = repmat(uA, nPaths, 1);

sqDt = sqrt(dt);

dW = randn(nStoch * nPaths, Nt);
dZ = 0.5 * sqDt^3 * (dW + 1/sqrt(3) * randn(nStoch * nPaths, Nt));
dW = sqDt * dW;

sig = repmat(sigma(:,1), nPaths, 1);
theta = repmat(sigma(:, 2), nPaths, 1);
sigthet = sig .* theta;
sigsq = sig .* sig;
sig3 = sigsq .* sig;
for n=1:Nt
    squn = sqrt(u(:,n));
    nZ = u(:,n) ~= 0;
    %% Implicit strong order 1 scheme
    u(:,n + 1) = u(:,n) + 0.5 * dt * theta .* (2 * uAtr - u(:, n)) ...
            + sig .* squn.*dW(:,n) + 0.25 * sigsq .* (dW(:,n).^2 - dt);
        
    %% Implicit strong order 1.5 scheme if not 0 in t(n)
    u(nZ,n + 1) = u(nZ,n + 1) + ...
        + (4* sigthet(nZ) .* (uAtr(nZ) - u(nZ, n)) - sig3(nZ))./(8 * squn(nZ)) ...
        .* (dW(nZ,n) * dt - dZ(nZ, n)) ...
        + sigsq(nZ) .* squn(nZ) .* (dZ(nZ, n) - 0.5 * dW(nZ, n) * dt);

    u(:,n + 1) = subplus(u(:,n + 1)) ./ (1 + 0.5 * dt * theta);
end
end

