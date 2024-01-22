%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'linear_sims';
M_.dynare_version = '5.4';
oo_.dynare_version = '5.4';
options_.dynare_version = '5.4';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(4,1);
M_.exo_names_tex = cell(4,1);
M_.exo_names_long = cell(4,1);
M_.exo_names(1) = {'eps_mu'};
M_.exo_names_tex(1) = {'eps\_mu'};
M_.exo_names_long(1) = {'eps_mu'};
M_.exo_names(2) = {'eps_mps'};
M_.exo_names_tex(2) = {'eps\_mps'};
M_.exo_names_long(2) = {'eps_mps'};
M_.exo_names(3) = {'eps_z'};
M_.exo_names_tex(3) = {'eps\_z'};
M_.exo_names_long(3) = {'eps_z'};
M_.exo_names(4) = {'eps_phi0'};
M_.exo_names_tex(4) = {'eps\_phi0'};
M_.exo_names_long(4) = {'eps_phi0'};
M_.endo_names = cell(63,1);
M_.endo_names_tex = cell(63,1);
M_.endo_names_long = cell(63,1);
M_.endo_names(1) = {'Q'};
M_.endo_names_tex(1) = {'Q'};
M_.endo_names_long(1) = {'Q'};
M_.endo_names(2) = {'lambda'};
M_.endo_names_tex(2) = {'lambda'};
M_.endo_names_long(2) = {'lambda'};
M_.endo_names(3) = {'R'};
M_.endo_names_tex(3) = {'R'};
M_.endo_names_long(3) = {'R'};
M_.endo_names(4) = {'X'};
M_.endo_names_tex(4) = {'X'};
M_.endo_names_long(4) = {'X'};
M_.endo_names(5) = {'r_u'};
M_.endo_names_tex(5) = {'r\_u'};
M_.endo_names_long(5) = {'r_u'};
M_.endo_names(6) = {'r_e'};
M_.endo_names_tex(6) = {'r\_e'};
M_.endo_names_long(6) = {'r_e'};
M_.endo_names(7) = {'teta'};
M_.endo_names_tex(7) = {'teta'};
M_.endo_names_long(7) = {'teta'};
M_.endo_names(8) = {'U'};
M_.endo_names_tex(8) = {'U'};
M_.endo_names_long(8) = {'U'};
M_.endo_names(9) = {'l_b'};
M_.endo_names_tex(9) = {'l\_b'};
M_.endo_names_long(9) = {'l_b'};
M_.endo_names(10) = {'l_g'};
M_.endo_names_tex(10) = {'l\_g'};
M_.endo_names_long(10) = {'l_g'};
M_.endo_names(11) = {'empl'};
M_.endo_names_tex(11) = {'empl'};
M_.endo_names_long(11) = {'empl'};
M_.endo_names(12) = {'mu'};
M_.endo_names_tex(12) = {'mu'};
M_.endo_names_long(12) = {'mu'};
M_.endo_names(13) = {'phi_0'};
M_.endo_names_tex(13) = {'phi\_0'};
M_.endo_names_long(13) = {'phi_0'};
M_.endo_names(14) = {'mps'};
M_.endo_names_tex(14) = {'mps'};
M_.endo_names_long(14) = {'mps'};
M_.endo_names(15) = {'z'};
M_.endo_names_tex(15) = {'z'};
M_.endo_names_long(15) = {'z'};
M_.endo_names(16) = {'EE'};
M_.endo_names_tex(16) = {'EE'};
M_.endo_names_long(16) = {'EE'};
M_.endo_names(17) = {'UE'};
M_.endo_names_tex(17) = {'UE'};
M_.endo_names_long(17) = {'UE'};
M_.endo_names(18) = {'AC'};
M_.endo_names_tex(18) = {'AC'};
M_.endo_names_long(18) = {'AC'};
M_.endo_names(19) = {'C'};
M_.endo_names_tex(19) = {'C'};
M_.endo_names_long(19) = {'C'};
M_.endo_names(20) = {'v'};
M_.endo_names_tex(20) = {'v'};
M_.endo_names_long(20) = {'v'};
M_.endo_names(21) = {'varphi'};
M_.endo_names_tex(21) = {'varphi'};
M_.endo_names_long(21) = {'varphi'};
M_.endo_names(22) = {'VU'};
M_.endo_names_tex(22) = {'VU'};
M_.endo_names_long(22) = {'VU'};
M_.endo_names(23) = {'Q_hat'};
M_.endo_names_tex(23) = {'Q\_hat'};
M_.endo_names_long(23) = {'Q_hat'};
M_.endo_names(24) = {'lambda_hat'};
M_.endo_names_tex(24) = {'lambda\_hat'};
M_.endo_names_long(24) = {'lambda_hat'};
M_.endo_names(25) = {'R_hat'};
M_.endo_names_tex(25) = {'R\_hat'};
M_.endo_names_long(25) = {'R_hat'};
M_.endo_names(26) = {'pi'};
M_.endo_names_tex(26) = {'pi'};
M_.endo_names_long(26) = {'pi'};
M_.endo_names(27) = {'W_hat'};
M_.endo_names_tex(27) = {'W\_hat'};
M_.endo_names_long(27) = {'W_hat'};
M_.endo_names(28) = {'X_hat'};
M_.endo_names_tex(28) = {'X\_hat'};
M_.endo_names_long(28) = {'X_hat'};
M_.endo_names(29) = {'r_u_hat'};
M_.endo_names_tex(29) = {'r\_u\_hat'};
M_.endo_names_long(29) = {'r_u_hat'};
M_.endo_names(30) = {'r_e_hat'};
M_.endo_names_tex(30) = {'r\_e\_hat'};
M_.endo_names_long(30) = {'r_e_hat'};
M_.endo_names(31) = {'mu_e_hat'};
M_.endo_names_tex(31) = {'mu\_e\_hat'};
M_.endo_names_long(31) = {'mu_e_hat'};
M_.endo_names(32) = {'teta_hat'};
M_.endo_names_tex(32) = {'teta\_hat'};
M_.endo_names_long(32) = {'teta_hat'};
M_.endo_names(33) = {'U_hat'};
M_.endo_names_tex(33) = {'U\_hat'};
M_.endo_names_long(33) = {'U_hat'};
M_.endo_names(34) = {'l_b_hat'};
M_.endo_names_tex(34) = {'l\_b\_hat'};
M_.endo_names_long(34) = {'l_b_hat'};
M_.endo_names(35) = {'l_g_hat'};
M_.endo_names_tex(35) = {'l\_g\_hat'};
M_.endo_names_long(35) = {'l_g_hat'};
M_.endo_names(36) = {'empl_hat'};
M_.endo_names_tex(36) = {'empl\_hat'};
M_.endo_names_long(36) = {'empl_hat'};
M_.endo_names(37) = {'mu_hat'};
M_.endo_names_tex(37) = {'mu\_hat'};
M_.endo_names_long(37) = {'mu_hat'};
M_.endo_names(38) = {'z_hat'};
M_.endo_names_tex(38) = {'z\_hat'};
M_.endo_names_long(38) = {'z_hat'};
M_.endo_names(39) = {'phi_0_hat'};
M_.endo_names_tex(39) = {'phi\_0\_hat'};
M_.endo_names_long(39) = {'phi_0_hat'};
M_.endo_names(40) = {'mps_hat'};
M_.endo_names_tex(40) = {'mps\_hat'};
M_.endo_names_long(40) = {'mps_hat'};
M_.endo_names(41) = {'INFL'};
M_.endo_names_tex(41) = {'INFL'};
M_.endo_names_long(41) = {'INFL'};
M_.endo_names(42) = {'EE_hat'};
M_.endo_names_tex(42) = {'EE\_hat'};
M_.endo_names_long(42) = {'EE_hat'};
M_.endo_names(43) = {'UE_hat'};
M_.endo_names_tex(43) = {'UE\_hat'};
M_.endo_names_long(43) = {'UE_hat'};
M_.endo_names(44) = {'AC_hat'};
M_.endo_names_tex(44) = {'AC\_hat'};
M_.endo_names_long(44) = {'AC_hat'};
M_.endo_names(45) = {'C_hat'};
M_.endo_names_tex(45) = {'C\_hat'};
M_.endo_names_long(45) = {'C_hat'};
M_.endo_names(46) = {'v_hat'};
M_.endo_names_tex(46) = {'v\_hat'};
M_.endo_names_long(46) = {'v_hat'};
M_.endo_names(47) = {'varphi_hat'};
M_.endo_names_tex(47) = {'varphi\_hat'};
M_.endo_names_long(47) = {'varphi_hat'};
M_.endo_names(48) = {'WplusLam'};
M_.endo_names_tex(48) = {'WplusLam'};
M_.endo_names_long(48) = {'WplusLam'};
M_.endo_names(49) = {'b_share'};
M_.endo_names_tex(49) = {'b\_share'};
M_.endo_names_long(49) = {'b_share'};
M_.endo_names(50) = {'H_hat'};
M_.endo_names_tex(50) = {'H\_hat'};
M_.endo_names_long(50) = {'H_hat'};
M_.endo_names(51) = {'TFP_hat'};
M_.endo_names_tex(51) = {'TFP\_hat'};
M_.endo_names_long(51) = {'TFP_hat'};
M_.endo_names(52) = {'EEf_hat'};
M_.endo_names_tex(52) = {'EEf\_hat'};
M_.endo_names_long(52) = {'EEf_hat'};
M_.endo_names(53) = {'UEf_hat'};
M_.endo_names_tex(53) = {'UEf\_hat'};
M_.endo_names_long(53) = {'UEf_hat'};
M_.endo_names(54) = {'AUX_ENDO_LAG_25_1'};
M_.endo_names_tex(54) = {'AUX\_ENDO\_LAG\_25\_1'};
M_.endo_names_long(54) = {'AUX_ENDO_LAG_25_1'};
M_.endo_names(55) = {'AUX_ENDO_LAG_25_2'};
M_.endo_names_tex(55) = {'AUX\_ENDO\_LAG\_25\_2'};
M_.endo_names_long(55) = {'AUX_ENDO_LAG_25_2'};
M_.endo_names(56) = {'AUX_ENDO_LAG_25_3'};
M_.endo_names_tex(56) = {'AUX\_ENDO\_LAG\_25\_3'};
M_.endo_names_long(56) = {'AUX_ENDO_LAG_25_3'};
M_.endo_names(57) = {'AUX_ENDO_LAG_25_4'};
M_.endo_names_tex(57) = {'AUX\_ENDO\_LAG\_25\_4'};
M_.endo_names_long(57) = {'AUX_ENDO_LAG_25_4'};
M_.endo_names(58) = {'AUX_ENDO_LAG_25_5'};
M_.endo_names_tex(58) = {'AUX\_ENDO\_LAG\_25\_5'};
M_.endo_names_long(58) = {'AUX_ENDO_LAG_25_5'};
M_.endo_names(59) = {'AUX_ENDO_LAG_25_6'};
M_.endo_names_tex(59) = {'AUX\_ENDO\_LAG\_25\_6'};
M_.endo_names_long(59) = {'AUX_ENDO_LAG_25_6'};
M_.endo_names(60) = {'AUX_ENDO_LAG_25_7'};
M_.endo_names_tex(60) = {'AUX\_ENDO\_LAG\_25\_7'};
M_.endo_names_long(60) = {'AUX_ENDO_LAG_25_7'};
M_.endo_names(61) = {'AUX_ENDO_LAG_25_8'};
M_.endo_names_tex(61) = {'AUX\_ENDO\_LAG\_25\_8'};
M_.endo_names_long(61) = {'AUX_ENDO_LAG_25_8'};
M_.endo_names(62) = {'AUX_ENDO_LAG_25_9'};
M_.endo_names_tex(62) = {'AUX\_ENDO\_LAG\_25\_9'};
M_.endo_names_long(62) = {'AUX_ENDO_LAG_25_9'};
M_.endo_names(63) = {'AUX_ENDO_LAG_25_10'};
M_.endo_names_tex(63) = {'AUX\_ENDO\_LAG\_25\_10'};
M_.endo_names_long(63) = {'AUX_ENDO_LAG_25_10'};
M_.endo_partitions = struct();
M_.param_names = cell(60,1);
M_.param_names_tex = cell(60,1);
M_.param_names_long = cell(60,1);
M_.param_names(1) = {'s'};
M_.param_names_tex(1) = {'s'};
M_.param_names_long(1) = {'s'};
M_.param_names(2) = {'y_b'};
M_.param_names_tex(2) = {'y\_b'};
M_.param_names_long(2) = {'y_b'};
M_.param_names(3) = {'s_g'};
M_.param_names_tex(3) = {'s\_g'};
M_.param_names_long(3) = {'s_g'};
M_.param_names(4) = {'xi_g'};
M_.param_names_tex(4) = {'xi\_g'};
M_.param_names_long(4) = {'xi_g'};
M_.param_names(5) = {'ups'};
M_.param_names_tex(5) = {'ups'};
M_.param_names_long(5) = {'ups'};
M_.param_names(6) = {'sigma'};
M_.param_names_tex(6) = {'sigma'};
M_.param_names_long(6) = {'sigma'};
M_.param_names(7) = {'iot'};
M_.param_names_tex(7) = {'iot'};
M_.param_names_long(7) = {'iot'};
M_.param_names(8) = {'varb'};
M_.param_names_tex(8) = {'varb'};
M_.param_names_long(8) = {'varb'};
M_.param_names(9) = {'psi'};
M_.param_names_tex(9) = {'psi'};
M_.param_names_long(9) = {'psi'};
M_.param_names(10) = {'xi'};
M_.param_names_tex(10) = {'xi'};
M_.param_names_long(10) = {'xi'};
M_.param_names(11) = {'zeta'};
M_.param_names_tex(11) = {'zeta'};
M_.param_names_long(11) = {'zeta'};
M_.param_names(12) = {'rho_mu'};
M_.param_names_tex(12) = {'rho\_mu'};
M_.param_names_long(12) = {'rho_mu'};
M_.param_names(13) = {'rho_mps'};
M_.param_names_tex(13) = {'rho\_mps'};
M_.param_names_long(13) = {'rho_mps'};
M_.param_names(14) = {'rho_z'};
M_.param_names_tex(14) = {'rho\_z'};
M_.param_names_long(14) = {'rho_z'};
M_.param_names(15) = {'rho_phi0'};
M_.param_names_tex(15) = {'rho\_phi0'};
M_.param_names_long(15) = {'rho_phi0'};
M_.param_names(16) = {'delta'};
M_.param_names_tex(16) = {'delta'};
M_.param_names_long(16) = {'delta'};
M_.param_names(17) = {'beta'};
M_.param_names_tex(17) = {'beta'};
M_.param_names_long(17) = {'beta'};
M_.param_names(18) = {'xi_b'};
M_.param_names_tex(18) = {'xi\_b'};
M_.param_names_long(18) = {'xi_b'};
M_.param_names(19) = {'y_g'};
M_.param_names_tex(19) = {'y\_g'};
M_.param_names_long(19) = {'y_g'};
M_.param_names(20) = {'b'};
M_.param_names_tex(20) = {'b'};
M_.param_names_long(20) = {'b'};
M_.param_names(21) = {'eta'};
M_.param_names_tex(21) = {'eta'};
M_.param_names_long(21) = {'eta'};
M_.param_names(22) = {'sz'};
M_.param_names_tex(22) = {'sz'};
M_.param_names_long(22) = {'sz'};
M_.param_names(23) = {'kap_s'};
M_.param_names_tex(23) = {'kap\_s'};
M_.param_names_long(23) = {'kap_s'};
M_.param_names(24) = {'c_p'};
M_.param_names_tex(24) = {'c\_p'};
M_.param_names_long(24) = {'c_p'};
M_.param_names(25) = {'phi_cap_ss'};
M_.param_names_tex(25) = {'phi\_cap\_ss'};
M_.param_names_long(25) = {'phi_cap_ss'};
M_.param_names(26) = {'rho_r'};
M_.param_names_tex(26) = {'rho\_r'};
M_.param_names_long(26) = {'rho_r'};
M_.param_names(27) = {'phi_pi'};
M_.param_names_tex(27) = {'phi\_pi'};
M_.param_names_long(27) = {'phi_pi'};
M_.param_names(28) = {'phi_u'};
M_.param_names_tex(28) = {'phi\_u'};
M_.param_names_long(28) = {'phi_u'};
M_.param_names(29) = {'phi_EE'};
M_.param_names_tex(29) = {'phi\_EE'};
M_.param_names_long(29) = {'phi_EE'};
M_.param_names(30) = {'phi_b_share'};
M_.param_names_tex(30) = {'phi\_b\_share'};
M_.param_names_long(30) = {'phi_b_share'};
M_.param_names(31) = {'mu_u'};
M_.param_names_tex(31) = {'mu\_u'};
M_.param_names_long(31) = {'mu_u'};
M_.param_names(32) = {'sigma_z'};
M_.param_names_tex(32) = {'sigma\_z'};
M_.param_names_long(32) = {'sigma_z'};
M_.param_names(33) = {'sigma_mps'};
M_.param_names_tex(33) = {'sigma\_mps'};
M_.param_names_long(33) = {'sigma_mps'};
M_.param_names(34) = {'sigma_mu'};
M_.param_names_tex(34) = {'sigma\_mu'};
M_.param_names_long(34) = {'sigma_mu'};
M_.param_names(35) = {'sigma_phi0'};
M_.param_names_tex(35) = {'sigma\_phi0'};
M_.param_names_long(35) = {'sigma_phi0'};
M_.param_names(36) = {'b_cap'};
M_.param_names_tex(36) = {'b\_cap'};
M_.param_names_long(36) = {'b_cap'};
M_.param_names(37) = {'varb_cap'};
M_.param_names_tex(37) = {'varb\_cap'};
M_.param_names_long(37) = {'varb_cap'};
M_.param_names(38) = {'prox_u'};
M_.param_names_tex(38) = {'prox\_u'};
M_.param_names_long(38) = {'prox_u'};
M_.param_names(39) = {'prox_e'};
M_.param_names_tex(39) = {'prox\_e'};
M_.param_names_long(39) = {'prox_e'};
M_.param_names(40) = {'Q_ss'};
M_.param_names_tex(40) = {'Q\_ss'};
M_.param_names_long(40) = {'Q_ss'};
M_.param_names(41) = {'lambda_ss'};
M_.param_names_tex(41) = {'lambda\_ss'};
M_.param_names_long(41) = {'lambda_ss'};
M_.param_names(42) = {'R_ss'};
M_.param_names_tex(42) = {'R\_ss'};
M_.param_names_long(42) = {'R_ss'};
M_.param_names(43) = {'X_ss'};
M_.param_names_tex(43) = {'X\_ss'};
M_.param_names_long(43) = {'X_ss'};
M_.param_names(44) = {'r_u_ss'};
M_.param_names_tex(44) = {'r\_u\_ss'};
M_.param_names_long(44) = {'r_u_ss'};
M_.param_names(45) = {'r_e_ss'};
M_.param_names_tex(45) = {'r\_e\_ss'};
M_.param_names_long(45) = {'r_e_ss'};
M_.param_names(46) = {'mu_e_ss'};
M_.param_names_tex(46) = {'mu\_e\_ss'};
M_.param_names_long(46) = {'mu_e_ss'};
M_.param_names(47) = {'teta_ss'};
M_.param_names_tex(47) = {'teta\_ss'};
M_.param_names_long(47) = {'teta_ss'};
M_.param_names(48) = {'U_ss'};
M_.param_names_tex(48) = {'U\_ss'};
M_.param_names_long(48) = {'U_ss'};
M_.param_names(49) = {'l_b_ss'};
M_.param_names_tex(49) = {'l\_b\_ss'};
M_.param_names_long(49) = {'l_b_ss'};
M_.param_names(50) = {'l_g_ss'};
M_.param_names_tex(50) = {'l\_g\_ss'};
M_.param_names_long(50) = {'l_g_ss'};
M_.param_names(51) = {'EE_ss'};
M_.param_names_tex(51) = {'EE\_ss'};
M_.param_names_long(51) = {'EE_ss'};
M_.param_names(52) = {'UE_ss'};
M_.param_names_tex(52) = {'UE\_ss'};
M_.param_names_long(52) = {'UE_ss'};
M_.param_names(53) = {'EU_ss'};
M_.param_names_tex(53) = {'EU\_ss'};
M_.param_names_long(53) = {'EU_ss'};
M_.param_names(54) = {'AC_ss'};
M_.param_names_tex(54) = {'AC\_ss'};
M_.param_names_long(54) = {'AC_ss'};
M_.param_names(55) = {'EEf_ss'};
M_.param_names_tex(55) = {'EEf\_ss'};
M_.param_names_long(55) = {'EEf_ss'};
M_.param_names(56) = {'b_share_ss'};
M_.param_names_tex(56) = {'b\_share\_ss'};
M_.param_names_long(56) = {'b_share_ss'};
M_.param_names(57) = {'v_ss'};
M_.param_names_tex(57) = {'v\_ss'};
M_.param_names_long(57) = {'v_ss'};
M_.param_names(58) = {'C_ss'};
M_.param_names_tex(58) = {'C\_ss'};
M_.param_names_long(58) = {'C_ss'};
M_.param_names(59) = {'varphi_ss'};
M_.param_names_tex(59) = {'varphi\_ss'};
M_.param_names_long(59) = {'varphi_ss'};
M_.param_names(60) = {'empl_ss'};
M_.param_names_tex(60) = {'empl\_ss'};
M_.param_names_long(60) = {'empl_ss'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 4;
M_.endo_nbr = 63;
M_.param_nbr = 60;
M_.orig_endo_nbr = 53;
M_.aux_vars(1).endo_index = 54;
M_.aux_vars(1).type = 1;
M_.aux_vars(1).orig_index = 26;
M_.aux_vars(1).orig_lead_lag = -1;
M_.aux_vars(1).orig_expr = 'pi(-1)';
M_.aux_vars(2).endo_index = 55;
M_.aux_vars(2).type = 1;
M_.aux_vars(2).orig_index = 26;
M_.aux_vars(2).orig_lead_lag = -2;
M_.aux_vars(2).orig_expr = 'AUX_ENDO_LAG_25_1(-1)';
M_.aux_vars(3).endo_index = 56;
M_.aux_vars(3).type = 1;
M_.aux_vars(3).orig_index = 26;
M_.aux_vars(3).orig_lead_lag = -3;
M_.aux_vars(3).orig_expr = 'AUX_ENDO_LAG_25_2(-1)';
M_.aux_vars(4).endo_index = 57;
M_.aux_vars(4).type = 1;
M_.aux_vars(4).orig_index = 26;
M_.aux_vars(4).orig_lead_lag = -4;
M_.aux_vars(4).orig_expr = 'AUX_ENDO_LAG_25_3(-1)';
M_.aux_vars(5).endo_index = 58;
M_.aux_vars(5).type = 1;
M_.aux_vars(5).orig_index = 26;
M_.aux_vars(5).orig_lead_lag = -5;
M_.aux_vars(5).orig_expr = 'AUX_ENDO_LAG_25_4(-1)';
M_.aux_vars(6).endo_index = 59;
M_.aux_vars(6).type = 1;
M_.aux_vars(6).orig_index = 26;
M_.aux_vars(6).orig_lead_lag = -6;
M_.aux_vars(6).orig_expr = 'AUX_ENDO_LAG_25_5(-1)';
M_.aux_vars(7).endo_index = 60;
M_.aux_vars(7).type = 1;
M_.aux_vars(7).orig_index = 26;
M_.aux_vars(7).orig_lead_lag = -7;
M_.aux_vars(7).orig_expr = 'AUX_ENDO_LAG_25_6(-1)';
M_.aux_vars(8).endo_index = 61;
M_.aux_vars(8).type = 1;
M_.aux_vars(8).orig_index = 26;
M_.aux_vars(8).orig_lead_lag = -8;
M_.aux_vars(8).orig_expr = 'AUX_ENDO_LAG_25_7(-1)';
M_.aux_vars(9).endo_index = 62;
M_.aux_vars(9).type = 1;
M_.aux_vars(9).orig_index = 26;
M_.aux_vars(9).orig_lead_lag = -9;
M_.aux_vars(9).orig_expr = 'AUX_ENDO_LAG_25_8(-1)';
M_.aux_vars(10).endo_index = 63;
M_.aux_vars(10).type = 1;
M_.aux_vars(10).orig_index = 26;
M_.aux_vars(10).orig_lead_lag = -10;
M_.aux_vars(10).orig_expr = 'AUX_ENDO_LAG_25_9(-1)';
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(4, 4);
M_.Correlation_matrix = eye(4, 4);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.orig_eq_nbr = 53;
M_.eq_nbr = 63;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 11;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 11;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 11;
M_.lead_lag_incidence = [
 0 20 0;
 0 21 0;
 0 22 0;
 0 23 0;
 0 24 0;
 0 25 0;
 0 26 0;
 0 27 0;
 0 28 0;
 0 29 0;
 0 30 0;
 0 31 0;
 0 32 0;
 0 33 0;
 0 34 0;
 0 35 0;
 0 36 0;
 0 37 0;
 0 38 0;
 0 39 0;
 0 40 0;
 0 41 0;
 0 42 0;
 0 43 83;
 1 44 0;
 2 45 84;
 0 46 85;
 0 47 0;
 0 48 0;
 0 49 0;
 3 50 0;
 0 51 0;
 4 52 0;
 5 53 0;
 0 54 0;
 0 55 0;
 6 56 0;
 7 57 0;
 8 58 0;
 9 59 0;
 0 60 0;
 0 61 0;
 0 62 0;
 0 63 0;
 0 64 0;
 0 65 0;
 0 66 0;
 0 67 0;
 0 68 0;
 0 69 0;
 0 70 0;
 0 71 0;
 0 72 0;
 10 73 0;
 11 74 0;
 12 75 0;
 13 76 0;
 14 77 0;
 15 78 0;
 16 79 0;
 17 80 0;
 18 81 0;
 19 82 0;]';
M_.nstatic = 42;
M_.nfwrd   = 2;
M_.npred   = 18;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 19;
M_.ndynamic   = 21;
M_.dynamic_tmp_nbr = [15; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
  3 , 'name' , '3' ;
  4 , 'name' , '4' ;
  5 , 'name' , 'mu_e_hat' ;
  6 , 'name' , '6' ;
  7 , 'name' , '7' ;
  8 , 'name' , '8' ;
  9 , 'name' , '9' ;
  10 , 'name' , '10' ;
  11 , 'name' , 'pi' ;
  12 , 'name' , '12' ;
  13 , 'name' , '13' ;
  14 , 'name' , '14' ;
  15 , 'name' , '15' ;
  16 , 'name' , '16' ;
  17 , 'name' , '17' ;
  18 , 'name' , '18' ;
  19 , 'name' , 'INFL' ;
  20 , 'name' , 'EE_hat' ;
  21 , 'name' , 'UE_hat' ;
  22 , 'name' , 'AC_hat' ;
  23 , 'name' , 'C_hat' ;
  24 , 'name' , 'v_hat' ;
  25 , 'name' , 'varphi_hat' ;
  26 , 'name' , 'WplusLam' ;
  27 , 'name' , 'b_share' ;
  28 , 'name' , 'H_hat' ;
  29 , 'name' , 'TFP_hat' ;
  30 , 'name' , 'EEf_hat' ;
  31 , 'name' , 'UEf_hat' ;
  32 , 'name' , 'Q' ;
  33 , 'name' , 'lambda' ;
  34 , 'name' , 'R' ;
  35 , 'name' , 'X' ;
  36 , 'name' , 'z' ;
  37 , 'name' , 'r_u' ;
  38 , 'name' , 'r_e' ;
  39 , 'name' , 'mu' ;
  40 , 'name' , 'teta' ;
  41 , 'name' , 'U' ;
  42 , 'name' , 'l_b' ;
  43 , 'name' , 'l_g' ;
  44 , 'name' , 'phi_0' ;
  45 , 'name' , 'mps' ;
  46 , 'name' , 'EE' ;
  47 , 'name' , 'UE' ;
  48 , 'name' , 'AC' ;
  49 , 'name' , 'C' ;
  50 , 'name' , 'v' ;
  51 , 'name' , 'varphi' ;
  52 , 'name' , 'empl' ;
  53 , 'name' , 'VU' ;
};
M_.mapping.Q.eqidx = [32 ];
M_.mapping.lambda.eqidx = [33 53 ];
M_.mapping.R.eqidx = [34 ];
M_.mapping.X.eqidx = [35 ];
M_.mapping.r_u.eqidx = [37 ];
M_.mapping.r_e.eqidx = [38 ];
M_.mapping.teta.eqidx = [40 ];
M_.mapping.U.eqidx = [41 ];
M_.mapping.l_b.eqidx = [42 ];
M_.mapping.l_g.eqidx = [43 ];
M_.mapping.empl.eqidx = [52 ];
M_.mapping.mu.eqidx = [39 ];
M_.mapping.phi_0.eqidx = [44 ];
M_.mapping.mps.eqidx = [45 ];
M_.mapping.z.eqidx = [36 ];
M_.mapping.EE.eqidx = [46 ];
M_.mapping.UE.eqidx = [47 ];
M_.mapping.AC.eqidx = [48 ];
M_.mapping.C.eqidx = [49 ];
M_.mapping.v.eqidx = [50 ];
M_.mapping.varphi.eqidx = [51 ];
M_.mapping.VU.eqidx = [53 ];
M_.mapping.Q_hat.eqidx = [6 11 13 29 32 ];
M_.mapping.lambda_hat.eqidx = [1 2 3 4 6 13 23 26 28 33 ];
M_.mapping.R_hat.eqidx = [1 14 34 ];
M_.mapping.pi.eqidx = [1 11 19 ];
M_.mapping.W_hat.eqidx = [2 3 4 26 ];
M_.mapping.X_hat.eqidx = [2 3 4 11 12 13 25 28 35 ];
M_.mapping.r_u_hat.eqidx = [3 7 9 12 20 21 37 ];
M_.mapping.r_e_hat.eqidx = [4 7 12 20 38 ];
M_.mapping.mu_e_hat.eqidx = [4 5 ];
M_.mapping.teta_hat.eqidx = [6 7 9 12 20 21 24 40 ];
M_.mapping.U_hat.eqidx = [6 7 8 9 10 12 13 14 20 24 31 41 ];
M_.mapping.l_b_hat.eqidx = [5 7 8 13 20 27 42 ];
M_.mapping.l_g_hat.eqidx = [8 43 ];
M_.mapping.empl_hat.eqidx = [5 10 27 29 30 52 ];
M_.mapping.mu_hat.eqidx = [6 13 16 23 39 ];
M_.mapping.z_hat.eqidx = [2 3 4 12 13 15 25 28 36 ];
M_.mapping.phi_0_hat.eqidx = [7 9 12 17 20 21 44 ];
M_.mapping.mps_hat.eqidx = [14 18 45 ];
M_.mapping.INFL.eqidx = [14 19 ];
M_.mapping.EE_hat.eqidx = [20 22 30 46 ];
M_.mapping.UE_hat.eqidx = [21 22 31 47 ];
M_.mapping.AC_hat.eqidx = [22 48 ];
M_.mapping.C_hat.eqidx = [23 49 ];
M_.mapping.v_hat.eqidx = [24 50 ];
M_.mapping.varphi_hat.eqidx = [25 51 ];
M_.mapping.WplusLam.eqidx = [26 ];
M_.mapping.b_share.eqidx = [14 27 ];
M_.mapping.H_hat.eqidx = [28 29 ];
M_.mapping.TFP_hat.eqidx = [29 ];
M_.mapping.EEf_hat.eqidx = [30 ];
M_.mapping.UEf_hat.eqidx = [31 ];
M_.mapping.eps_mu.eqidx = [16 ];
M_.mapping.eps_mps.eqidx = [18 ];
M_.mapping.eps_z.eqidx = [15 ];
M_.mapping.eps_phi0.eqidx = [17 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [25 26 31 33 34 37 38 39 40 54 55 56 57 58 59 60 61 62 63 ];
M_.exo_names_orig_ord = [1:4];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(63, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(4, 1);
M_.params = NaN(60, 1);
M_.endo_trends = struct('deflator', cell(63, 1), 'log_deflator', cell(63, 1), 'growth_factor', cell(63, 1), 'log_growth_factor', cell(63, 1));
M_.NNZDerivatives = [204; -1; -1; ];
M_.static_tmp_nbr = [15; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 0.99;
s = M_.params(1);
M_.params(2) = 1;
y_b = M_.params(2);
sg = 0.01; 
M_.params(4) = 0.5;
xi_g = M_.params(4);
M_.params(5) = 0.67;
ups = M_.params(5);
M_.params(6) = 0.33;
sigma = M_.params(6);
M_.params(7) = 2.35;
iot = M_.params(7);
M_.params(8) = 0.37;
varb = M_.params(8);
M_.params(9) = 0.38;
psi = M_.params(9);
M_.params(10) = 0.37;
xi = M_.params(10);
M_.params(11) = 0.11;
zeta = M_.params(11);
M_.params(51) = 0.024;
EE_ss = M_.params(51);
M_.params(52) = 0.25;
UE_ss = M_.params(52);
M_.params(53) = 0.01314;
EU_ss = M_.params(53);
M_.params(54) = M_.params(51)/M_.params(52);
AC_ss = M_.params(54);
M_.params(48) = M_.params(53)/(M_.params(52)+M_.params(53));
U_ss = M_.params(48);
M_.params(60) = 1-M_.params(48);
empl_ss = M_.params(60);
M_.params(55) = M_.params(51)*M_.params(60);
EEf_ss = M_.params(55);
M_.params(47) = 1;
teta_ss = M_.params(47);
M_.params(25) = 1;
phi_cap_ss = M_.params(25);
M_.params(21) = 6;
eta = M_.params(21);
M_.params(17) = 0.995;
beta = M_.params(17);
M_.params(12) = 0.94;
rho_mu = M_.params(12);
M_.params(14) = 0.494;
rho_z = M_.params(14);
M_.params(15) = 0.5;
rho_phi0 = M_.params(15);
M_.params(42) = 1/M_.params(17)-1;
R_ss = M_.params(42);
M_.params(32) = 0.00134689;
sigma_z = M_.params(32);
M_.params(33) = 0.000529;
sigma_mps = M_.params(33);
M_.params(34) = 0.000529;
sigma_mu = M_.params(34);
M_.params(35) = 0.000529;
sigma_phi0 = M_.params(35);
M_.params(13) = 0.93;
rho_mps = M_.params(13);
M_.params(27) = 1.17;
phi_pi = M_.params(27);
M_.params(28) = (-0.05);
phi_u = M_.params(28);
M_.params(29) = 0.05;
phi_EE = M_.params(29);
M_.params(30) = 0.05;
phi_b_share = M_.params(30);
M_.params(26) = 0.85;
rho_r = M_.params(26);
M_.params(19) = M_.params(2)*(1+sg);
y_g = M_.params(19);
M_.params(18) = 1-M_.params(4);
xi_b = M_.params(18);
M_.params(16) = M_.params(4)*(M_.params(51)+M_.params(53))/(M_.params(4)+M_.params(4)*M_.params(18)-M_.params(54));
delta = M_.params(16);
M_.params(22) = (1-M_.params(53)/M_.params(16))/M_.params(52);
sz = M_.params(22);
r_ratio = (s*AC_ss*(1-delta)*UE_ss)/(EE_ss-delta+EU_ss); 
M_.params(49) = M_.params(18)*(1-M_.params(48))*M_.params(16)/(M_.params(16)+M_.params(4)*(M_.params(51)-M_.params(52)*M_.params(16)*M_.params(22))/M_.params(54));
l_b_ss = M_.params(49);
M_.params(50) = 1-M_.params(48)-M_.params(49);
l_g_ss = M_.params(50);
M_.params(56) = M_.params(49)/M_.params(50);
b_share_ss = M_.params(56);
M_.params(31) = M_.params(4)*M_.params(19)^(1+M_.params(10))+M_.params(18)*M_.params(2)^(1+M_.params(10));
mu_u = M_.params(31);
M_.params(46) = M_.params(4)*M_.params(49)/(1-M_.params(48))*(M_.params(19)^(1+M_.params(10))-M_.params(2)^(1+M_.params(10)));
mu_e_ss = M_.params(46);
M_.params(41) = lambda_solver(M_.params(5),M_.params(6),M_.params(21),M_.params(17),M_.params(7),M_.params(10),M_.params(8),M_.params(1),M_.params(49),M_.params(50),M_.params(46),M_.params(2),M_.params(19),M_.params(16),M_.params(22),M_.params(51),M_.params(52),M_.params(53),M_.params(54),M_.params(48));
lambda_ss = M_.params(41);
M_.params(40) = Q_solver(M_.params(5),M_.params(6),M_.params(21),M_.params(17),M_.params(7),M_.params(10),M_.params(8),M_.params(1),M_.params(49),M_.params(50),M_.params(46),M_.params(2),M_.params(19),M_.params(16),M_.params(22),M_.params(51),M_.params(52),M_.params(53),M_.params(54),M_.params(48));
Q_ss = M_.params(40);
M_.params(43) = M_.params(5)*(M_.params(21)-1)/M_.params(21)*M_.params(40)^((M_.params(5)-1)/M_.params(5));
X_ss = M_.params(43);
M_.params(37) = 1/(1+M_.params(10))*(M_.params(43)*M_.params(41)/M_.params(8))^M_.params(10);
varb_cap = M_.params(37);
M_.params(36) = M_.params(37)*M_.params(31)-M_.params(46)*M_.params(37)*(M_.params(52)*M_.params(54)*M_.params(1)*(1-M_.params(16))/(M_.params(53)+M_.params(51)-M_.params(16)))^M_.params(7);
b_cap = M_.params(36);
M_.params(20) = M_.params(41)*M_.params(43)*M_.params(36);
b = M_.params(20);
M_.params(38) = M_.params(17)/(1-M_.params(17)*(1-M_.params(16)))*(M_.params(37)*M_.params(31)-M_.params(36));
prox_u = M_.params(38);
M_.params(39) = M_.params(46)*M_.params(37)*M_.params(17)/(1-M_.params(17)*(1-M_.params(16)));
prox_e = M_.params(39);
M_.params(24) = (M_.params(40)-M_.params(41)^(-M_.params(6)))/(M_.params(48)+(1-M_.params(48))*(M_.params(16)*M_.params(22)+M_.params(1)*(1-M_.params(16))));
c_p = M_.params(24);
M_.params(44) = M_.params(52);
r_u_ss = M_.params(44);
M_.params(45) = (M_.params(51)-M_.params(52)*M_.params(16)*M_.params(22))/(M_.params(1)*M_.params(54)*(1-M_.params(16)));
r_e_ss = M_.params(45);
M_.params(23) = M_.params(38)/M_.params(44)^M_.params(7);
kap_s = M_.params(23);
M_.params(58) = M_.params(41)^(-M_.params(6));
C_ss = M_.params(58);
M_.params(57) = M_.params(48)+(1-M_.params(48))*(M_.params(16)*M_.params(22)+M_.params(1)*(1-M_.params(16)));
v_ss = M_.params(57);
M_.params(59) = M_.params(43);
varphi_ss = M_.params(59);
VU_ss = (beta*b)/(lambda_ss*(1-beta*(1-delta)));
steady;
resid;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(34);
M_.Sigma_e(2, 2) = M_.params(33);
M_.Sigma_e(3, 3) = M_.params(32);
options_.impulse_responses.plot_threshold = 0;
options_.irf = 0;
options_.order = 1;
options_.periods = 10000;
var_list_ = {'mps';'EEf_hat';'EE_hat';'U_hat';'C_hat';'b_share'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'linear_sims_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'linear_sims_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'linear_sims_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'linear_sims_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'linear_sims_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'linear_sims_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'linear_sims_results.mat'], 'oo_recursive_', '-append');
end
disp('Note: 31 warning(s) encountered in the preprocessor')
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end