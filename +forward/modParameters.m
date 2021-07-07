classdef modParameters
    % Class containing the model Parameters
    %
    % Calling modParameters() will create an instance with default
    % parameters
    properties
        
        Dose1_2_5(1,1) double {mustBeNonnegative, mustBeFinite} = 2.5000
        K_mran(1,1) double {mustBeNonnegative, mustBeFinite} = 223.5815
        R1_total(1,1) double {mustBeNonnegative, mustBeFinite} = 34.8732
        R2_total(1,1) double {mustBeNonnegative, mustBeFinite} =34.9992
        S2_export_from_nuc(1,1) double {mustBeNonnegative, mustBeFinite} = 1.9990
        S2_import_to_nuc(1,1) double {mustBeNonnegative, mustBeFinite} = 0.3717
        S2_total(1,1) double {mustBeNonnegative, mustBeFinite} = 875.9947
        S4_export_from_nuc(1,1) double {mustBeNonnegative, mustBeFinite} = 0.1969
        S4_import_to_nuc(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0371
        S4_total(1,1) double {mustBeNonnegative, mustBeFinite} = 629.9245
        Trimer_import_to_nuc(1,1) double {mustBeNonnegative, mustBeFinite} = 0.1440
        export_cytoplasm(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0304
        hill_fact1(1,1) double {mustBeNonnegative, mustBeFinite} = 3.9865
        index_active_Rec_internalize(1,1) double {mustBeNonnegative, mustBeFinite} = 0.5376
        index_induced_R2_deg(1,1) double {mustBeNonnegative, mustBeFinite} = 1.0004
        index_induced_ligand_deg(1,1) double {mustBeNonnegative, mustBeFinite} = 2.7211
        index_k_out_1_relative_speed(1,1) double {mustBeNonnegative, mustBeFinite} = 0.5668
        index_k_out_2_relative_speed(1,1) double {mustBeNonnegative, mustBeFinite} = 0.2457
        index_kb_R1(1,1) double {mustBeNonnegative, mustBeFinite} = 1.6156
        index_kb_R2(1,1) double {mustBeNonnegative, mustBeFinite} = 8.7047
        index_kb_homotrimer(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0164
        index_kf_homotrimer(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(0.6519)
        index_seq_kb(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-3.0563)
        index_trimer_dephos(1,1) double {mustBeNonnegative, mustBeFinite} = 2.6103
        k1(1,1) double {mustBeNonnegative, mustBeFinite} = 0.4700
        k2(1,1) double {mustBeNonnegative, mustBeFinite} = 1.1750
        k3(1,1) double {mustBeNonnegative, mustBeFinite} = 2.3500
        k4(1,1) double {mustBeNonnegative, mustBeFinite} = 11.7500
        k5(1,1) double {mustBeNonnegative, mustBeFinite} = 47.0000
        kS2_effective_1(1,1) double {mustBeNonnegative, mustBeFinite} = 0.8500
        kS2_effective_2(1,1) double {mustBeNonnegative, mustBeFinite} = 0.9998
        kS4_effective_1(1,1) double {mustBeNonnegative, mustBeFinite} = 1.0000
        kS4_effective_2(1,1) double {mustBeNonnegative, mustBeFinite} = 1.0000
        k_Dephos(1,1) double {mustBeNonnegative, mustBeFinite} = 0.1394
        k_S7_protein(1,1) double {mustBeNonnegative, mustBeFinite} = 0.2183
        k_disso_Active_Rec(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0010
        k_in_1(1,1) double {mustBeNonnegative, mustBeFinite} = 0.2007
        k_in_2(1,1) double {mustBeNonnegative, mustBeFinite} = 0.6257
        k_induced_S7_production(1,1) double {mustBeNonnegative, mustBeFinite} = 4.9997
        k_medium(1,1) double {mustBeNonnegative, mustBeFinite} = 2018
        k_phosphorylation(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0702
        kb_trimmer(1,1) double {mustBeNonnegative, mustBeFinite} = 1.6078
        kdeg_R1(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0588
        kdeg_R2(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0010
        kdeg_S2(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0005
        kdeg_S4(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0048
        kdeg_S7(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0010
        kf_R1_activation(1,1) double {mustBeNonnegative, mustBeFinite} = 4.9855
        kf_R2_activation(1,1) double {mustBeNonnegative, mustBeFinite} = 0.1*4.9556
        kf_Seq_S7_Rec(1,1) double {mustBeNonnegative, mustBeFinite} = 0.8825
        kf_trimmer(1,1) double {mustBeNonnegative, mustBeFinite} = 4.8142
        kin_deg_Ligand(1,1) double {mustBeNonnegative, mustBeFinite} = 0.7201
        kmRNA1deg_S7(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0053
        kmRNAdeg_S7(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0990
        mRNA_prod(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0131
        offset_DRB_100_old(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-0.5432)
        offset_S2_Rec_Ini(1,1) double {mustBeNonnegative, mustBeFinite} = 0.0395
        offset_S4_Rec_Ini(1,1) double {mustBeNonnegative, mustBeFinite} = 0.2489
        offset_S7_mRNA(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-0.0262)
        offset_restimulation(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-0.4623)
        scale_DRB_100_old(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-0.0605)
        scale_S2_Rec_Ini(1,1) double {mustBeNonnegative, mustBeFinite} = 0.9471
        scale_S4_Rec_Ini(1,1) double {mustBeNonnegative, mustBeFinite} = 0.6210
        scale_S7_mRNA(1,1) double {mustBeNonnegative, mustBeFinite} = 0.1200
        scale_restimulation(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-0.3583)
        sd_Rec1(1,1) double {mustBeNonnegative, mustBeFinite} = 0.4000
        sd_Rec2(1,1) double {mustBeNonnegative, mustBeFinite} = 0.3000
        sd_S2_R_INI(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-1.6130)
        sd_S4_R_INI(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-1.3162)
        sd_S4_nExpID17(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-1.5068)
        sd_S4_nExpID7(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-1.0000)
        sd_S7_mRNA(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-0.0255)
        sd_pS2_nExpID12(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-1.2977)
        sd_pS2_nExpID30(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-1.5737)
        sd_pS2_nExpID31(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-1.6533)
        sd_pS2_nExpID32(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-1.5889)
        sd_pS2_nExpID37(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-0.9225)
        sd_pS2_nExpID7(1,1) double {mustBeNonnegative, mustBeFinite} = 10^(-1.0000)
        time1_2_5(1,1) double {mustBeNonnegative, mustBeFinite} = 10.0000
        
        alterPassiveDeg = []
        
        
        preDistribute = []
        
        
        test = false
    end
    methods(Static)
        function pfields = getFieldNames()
            % Output the parameter names
            
            import forward.modParameters % strange but necessary
            pfields = cell(102, 1);
            fnames =fieldnames(modParameters);
            pfields(24:102) = fnames(1:79);
        end
        
    end
    methods
        function p = unwrap(P)
            % Outputs parameter values in default order
            
            p = [ P.Dose1_2_5; P.K_mran; P.R1_total; P.R2_total; P.S2_export_from_nuc; ...
                P.S2_import_to_nuc; P.S2_total; P.S4_export_from_nuc; P.S4_import_to_nuc; P.S4_total; ...
                P.Trimer_import_to_nuc; P.export_cytoplasm; P.hill_fact1; P.index_active_Rec_internalize; ...
                P.index_induced_R2_deg; P.index_induced_ligand_deg; P.index_k_out_1_relative_speed; ...
                P.index_k_out_2_relative_speed; P.index_kb_R1; P.index_kb_R2; P.index_kb_homotrimer; ...
                P.index_kf_homotrimer; P.index_seq_kb; P.index_trimer_dephos; P.k1; P.k2; P.k3; P.k4; P.k5; ...
                P.kS2_effective_1; P.kS2_effective_2; P.kS4_effective_1; P.kS4_effective_2; P.k_Dephos; ...
                P.k_S7_protein; P.k_disso_Active_Rec; P.k_in_1; P.k_in_2; P.k_induced_S7_production; ...
                P.k_medium; P.k_phosphorylation; P.kb_trimmer; P.kdeg_R1; P.kdeg_R2; P.kdeg_S2; P.kdeg_S4; ...
                P.kdeg_S7; P.kf_R1_activation; P.kf_R2_activation; P.kf_Seq_S7_Rec; P.kf_trimmer; P.kin_deg_Ligand; ...
                P.kmRNA1deg_S7; P.kmRNAdeg_S7; P.mRNA_prod; P.offset_DRB_100_old; P.offset_S2_Rec_Ini; P.offset_S4_Rec_Ini; ...
                P.offset_S7_mRNA; P.offset_restimulation; P.scale_DRB_100_old; P.scale_S2_Rec_Ini; P.scale_S4_Rec_Ini; ...
                P.scale_S7_mRNA; P.scale_restimulation; P.sd_Rec1; P.sd_Rec2; P.sd_S2_R_INI; P.sd_S4_R_INI; P.sd_S4_nExpID17; ...
                P.sd_S4_nExpID7; P.sd_S7_mRNA; P.sd_pS2_nExpID12; P.sd_pS2_nExpID30; P.sd_pS2_nExpID31; P.sd_pS2_nExpID32; ...
                P.sd_pS2_nExpID37; P.sd_pS2_nExpID7; P.time1_2_5];
        end
    end
    
end