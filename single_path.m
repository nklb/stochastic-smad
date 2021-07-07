% Simulation of a single cell
close all
import forward.modParameters forward.getInitialData ...
    forward.compPathsVectorized

%% setup deterministic and stochastic parameters
P = modParameters();
P.preDistribute = 120;
sigma = zeros(102,2);

%% select from the five stochastic blockmodels degradation, endosomal
%% internalization, receptor-ligand and synthesis
model = 'internalization';

switch model
    case 'degradation'
        stoch_comp = [66; 67; 70]; 
        sigma(stoch_comp,1) = [1.58e-6;   2.04e-1;   1.92e-2]; 
        sigma(stoch_comp,2) = [ 6.76e-5;   9.23e-2;   7.90e-4];
        P.kdeg_R1 = 3.13e-3; P.kdeg_R2 = 1.01e-6; P.kdeg_S7 = 1.77e-3;    
    case 'endosomal'
        stoch_comp =  59; 
        sigma(stoch_comp,1) = [2.50e-1];
        sigma(stoch_comp,2) = [1.49e-2];
        P.k_disso_Active_Rec = 1.00e-2; 
    case 'internalization'
        stoch_comp = [37; 60; 61];
        sigma(stoch_comp,1) = [3.14e-2;   4.70e-2;   8.37e-2];
        sigma(stoch_comp,2) = [7.34e-3;   4.07e-3;   1.60e-2];
        P.index_active_Rec_internalize = 8.12e-1; P.k_in_1 = 3.78e-1; P.k_in_2 = 3.55e-1;
    case 'receptor-ligand'
        stoch_comp =  [71; 72];
        sigma(stoch_comp,1) = [1.13e-1;   1.54e-2];
        sigma(stoch_comp,2) = [3.33e-3;   5.37e-7];
        P.kf_R1_activation = 4.874; P.kf_R2_activation = 5.885;
    case 'synthesis'
        stoch_comp =  [26; 27; 78];
        sigma(stoch_comp,1) = [1.61e-8;   7.52e-7;   6.28e-5];
        sigma(stoch_comp,2) = [7.00e-15;   9.10e-11;   9.37e-6];
        P.R1_total = 27.364; P.R2_total = 33.614; P.mRNA_prod = 1.21e-3;
    otherwise
        error('The stochastic model %s is not defined!', model);
end

%% set options
% dose of TGFb (0: 0 pM, 1: 1 pM, 2: 2.5 pM,  3: 5 pM, 4: 25 pM, 
%               5: 100 pM, 6: 2 x 5 pM)
dose = 5; 

% remove negative paths
remNeg = true;

% fix seed (or not)
sameSeed = false;

% compute objective function
obj_function = true;

% add control paths
controls = true;

% time increment (mins)
dt = 0.6;

% final time (mins)
T = 1440;

%% simulation
[compPaths] = compPathsVectorized(1, P, sigma, dose, remNeg, sameSeed, ...
    dt, T);

%% visualisation
figure()
[t, smad2] = compPaths.getPrepData(controls);
plot(t/60, smad2);
xlim([0, 24]);
xlabel('time (h)');
ylabel('nuc/cyt SMAD2')