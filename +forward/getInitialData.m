function y0 = getInitialData(P)
% Provides initial data of concentrations and parameters for a given set of
% model parameters
%
% Pars:
%   P (obj:modParameters):      deterministic parameters used for the
%                               simulations. 
%
% Returns:
% y0 (vector[102]):             initial value for the time integration 
%                               including concentrations in the entries 
%                               1:23 and parameters in the entries 24:102

y0 = zeros(102,1);

%% define initial concentrations
y0(12) = (((P.S2_total * ((P.S2_import_to_nuc*P.kdeg_S2 + P.S2_export_from_nuc*P.kdeg_S2 + P.kdeg_S2^2)/((P.S2_export_from_nuc + P.kdeg_S2) + 2*P.S2_import_to_nuc)))*(P.S2_export_from_nuc + P.kdeg_S2))/(P.S2_import_to_nuc*P.kdeg_S2 + P.S2_export_from_nuc*P.kdeg_S2 + P.kdeg_S2^2));
y0(22) = ((2*P.S2_import_to_nuc*(P.S2_total * ((P.S2_import_to_nuc*P.kdeg_S2 + P.S2_export_from_nuc*P.kdeg_S2 + P.kdeg_S2^2)/((P.S2_export_from_nuc + P.kdeg_S2) + 2*P.S2_import_to_nuc))))/(P.S2_import_to_nuc*P.kdeg_S2 + P.S2_export_from_nuc*P.kdeg_S2 + P.kdeg_S2^2));
y0(13) = (((P.S4_total * ((P.S4_import_to_nuc*P.kdeg_S4 + P.S4_export_from_nuc*P.kdeg_S4 + P.kdeg_S4^2)/((P.S4_export_from_nuc + P.kdeg_S4) + 2*P.S4_import_to_nuc)))*(P.S4_export_from_nuc + P.kdeg_S4))/(P.S4_import_to_nuc*P.kdeg_S4 + P.S4_export_from_nuc*P.kdeg_S4 + P.kdeg_S4^2));
y0(23) = ((2*P.S4_import_to_nuc*(P.S4_total * ((P.S4_import_to_nuc*P.kdeg_S4 + P.S4_export_from_nuc*P.kdeg_S4 + P.kdeg_S4^2)/((P.S4_export_from_nuc + P.kdeg_S4) + 2*P.S4_import_to_nuc))))/(P.S4_import_to_nuc*P.kdeg_S4 + P.S4_export_from_nuc*P.kdeg_S4 + P.kdeg_S4^2));
y0(18) = ((P.export_cytoplasm*P.k_S7_protein*P.mRNA_prod)/(2*P.kdeg_S7*P.kmRNAdeg_S7*(P.export_cytoplasm + P.kmRNA1deg_S7)));
y0(16) = ((P.export_cytoplasm*P.mRNA_prod)/(2*P.kmRNAdeg_S7*(P.export_cytoplasm + P.kmRNA1deg_S7)));
y0(17) = ( P.mRNA_prod/(P.export_cytoplasm + P.kmRNA1deg_S7));
y0(5) = 0.00000000000001;
y0(3) = ((P.R1_total * ( P.k_in_1*(P.kdeg_R1)^2 )/((P.index_k_out_1_relative_speed * P.k_in_1) * P.kdeg_R1 + P.kdeg_R1^2 + P.k_in_1*P.kdeg_R1))/P.kdeg_R1);
y0(1) = (((P.R1_total * ( P.k_in_1*(P.kdeg_R1)^2 )/((P.index_k_out_1_relative_speed * P.k_in_1) * P.kdeg_R1 + P.kdeg_R1^2 + P.k_in_1*P.kdeg_R1))*((P.index_k_out_1_relative_speed * P.k_in_1) + P.kdeg_R1))/(P.k_in_1*P.kdeg_R1));
y0(4) = ((P.R2_total * ( P.k_in_2*(P.kdeg_R2)^2 )/((P.index_k_out_2_relative_speed * P.k_in_2) * P.kdeg_R2 + P.kdeg_R2^2 + P.k_in_2*P.kdeg_R2))/P.kdeg_R2);
y0(2) = (((P.R2_total * ( P.k_in_2*(P.kdeg_R2)^2 )/((P.index_k_out_2_relative_speed * P.k_in_2) * P.kdeg_R2 + P.kdeg_R2^2 + P.k_in_2*P.kdeg_R2))*((P.index_k_out_2_relative_speed * P.k_in_2) + P.kdeg_R2))/(P.k_in_2*P.kdeg_R2));


y0(24:end) = P.unwrap();
