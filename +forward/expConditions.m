classdef expConditions
    % Class containing the experimental conditions
    %
    % These conditions depend on both the deterministic parameters and the
    % experiment case.
    properties
        Reaction_kickstart(1,1) function_handle = @(t)0
        BolusInjection(1,1) function_handle = @(t)0
        passiveDeg(1,1) double {mustBeNonnegative, mustBeFinite}
        
        K_BolusInjection1(1,1) double {mustBeNonnegative, mustBeFinite} = 0
        K_BolusInjection2(1,1) double {mustBeNonnegative, mustBeFinite} = 0
        K_BolusInjection3(1,1) double {mustBeNonnegative, mustBeFinite} = 0
        K_DRB_index(1,1) double {mustBeNonnegative, mustBeFinite} = 1
        K_Rec_inhibitor(1,1) double {mustBeNonnegative, mustBeFinite} = 0
        K_S7_KD_index(1,1) double {mustBeNonnegative, mustBeFinite} = 1
        K_Wash(1,1) double {mustBeNonnegative, mustBeFinite} = 0
        K_Wash2(1,1) double {mustBeNonnegative, mustBeFinite} = 0
        K_cyclo_index(1,1) double {mustBeNonnegative, mustBeFinite} = 1
        K_effective(1,1) double {mustBeNonnegative, mustBeFinite} = 0
        K_induced_R1_production(1,1) double {mustBeNonnegative, mustBeFinite} = 0
        K_induced_R2_production(1,1) double {mustBeNonnegative, mustBeFinite} = 0
        K_k_R1_production(1,1) double {mustBeNonnegative, mustBeFinite}
        K_k_R2_production(1,1) double {mustBeNonnegative, mustBeFinite}
        K_k_S2_production(1,1) double {mustBeNonnegative, mustBeFinite}
        K_k_S4_production(1,1) double {mustBeNonnegative, mustBeFinite}
        K_kb_R1_activation(1,1) double {mustBeNonnegative, mustBeFinite}
        K_kb_R2_activation(1,1) double {mustBeNonnegative, mustBeFinite}
        K_kb_Seq_S7_Rec(1,1) double {mustBeNonnegative, mustBeFinite}
        
        TGFin = []
        TGFinT = []
        
        k_inj1 = 0
        k_inj2 = 0
    end
    methods
        function E = expConditions(P, dose)
            % Constructor function
            %
            % Args:
            %   P (obj:modParameters): Deterministic parameters used for the
            %    simulations
            %   dose(integer between 0 and 6): Defines the experiment case. 
            import burstDetection.getRootDir
            
            if ~isa(P, 'forward.modParameters')
                error('Input P must be an object of the class modParameters.')
            end
            E.Reaction_kickstart = @(t)1000./(1000 + (exp(-(t - 0.1)).^100));
            
            kond7 = 1;
            kond8 = 1;
            
            inj1 = {0, P.k1, P.k2, P.k3, P.k4, P.k5, P.k2, kond7*P.k2,kond8*P.k2};
            inj2 = {0, 0, 0, 0, 0, 0, P.k2,kond7*P.k2,kond8*P.k2};
            
            passiveDeg = {10.665, 10.665, 10.355, 10.355, 3.16, ...
                0.586, 10.664,0.665,10.664,0.665};
            E.k_inj1 = inj1{dose + 1};
            E.k_inj2 = inj2{dose + 1};
            E.BolusInjection = @(t, k_inj1, k_inj2) k_inj2 *  (1000./(1000 + (exp(-(t - 360.1)).^100)) .* 1000./(1000 + exp(t- 362.1).^100)) ...
                + k_inj1 *  (1000./(1000 + (exp(-(t - 0.1)).^100)) .* 1000./(1000 + exp(t- 2.1).^100));
            E.passiveDeg = passiveDeg{dose + 1};
            
            %% set experimental parameters
            
            E.K_k_R1_production = (P.R1_total*(P.k_in_1*(P.kdeg_R1)^2)/((P.index_k_out_1_relative_speed*P.k_in_1)*P.kdeg_R1+P.kdeg_R1^2+P.k_in_1*P.kdeg_R1));
            E.K_k_R2_production = (P.R2_total*(P.k_in_2*(P.kdeg_R2)^2)/((P.index_k_out_2_relative_speed*P.k_in_2)*P.kdeg_R2+P.kdeg_R2^2+P.k_in_2*P.kdeg_R2));
            E.K_k_S2_production = (P.S2_total*((P.S2_import_to_nuc*P.kdeg_S2+P.S2_export_from_nuc*P.kdeg_S2+P.kdeg_S2^2)/((P.S2_export_from_nuc+P.kdeg_S2)+2*P.S2_import_to_nuc)));
            E.K_k_S4_production = (P.S4_total*((P.S4_import_to_nuc*P.kdeg_S4+P.S4_export_from_nuc*P.kdeg_S4+P.kdeg_S4^2)/((P.S4_export_from_nuc+P.kdeg_S4)+2*P.S4_import_to_nuc)));
            
            E.K_kb_R1_activation = (P.kf_R1_activation * P.index_kb_R1);
            E.K_kb_R2_activation = (P.kf_R2_activation * P.index_kb_R2);
            E.K_kb_Seq_S7_Rec = (P.kf_Seq_S7_Rec * P.index_seq_kb);
            
            if ~isempty(P.alterPassiveDeg)
                E.passiveDeg = P.alterPassiveDeg;
            end
            

        end
    end
end
