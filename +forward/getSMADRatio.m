function SMAD_ratio = getSMADRatio(Y)
% Derives nuc/cyt SMAD2 ratio from dynamic variables
%
% Pars:
%   Y (102xN array):            dynamic concentrations at N time instances
%
% Returns:
% SMAD_ratio (Nx1 array):       nuc/cyt SMAD2 ratio at N time instances

SMAD_ratio = (Y(19, :) + 2*Y(20, :) + 3*Y(21, :) + Y(22, :)) ...
    ./(Y(11, :) + Y(12, :) + 2*Y(14, :) + 3*Y(15,:));
