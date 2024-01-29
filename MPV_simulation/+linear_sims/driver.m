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
M_.exo_names = cell(3,1);
M_.exo_names_tex = cell(3,1);
M_.exo_names_long = cell(3,1);
M_.exo_names(1) = {'eps_mu'};
M_.exo_names_tex(1) = {'eps\_mu'};
M_.exo_names_long(1) = {'eps_mu'};
M_.exo_names(2) = {'eps_z'};
M_.exo_names_tex(2) = {'eps\_z'};
M_.exo_names_long(2) = {'eps_z'};
M_.exo_names(3) = {'eps_mps'};
M_.exo_names_tex(3) = {'eps\_mps'};
M_.exo_names_long(3) = {'eps_mps'};
M_.endo_names = cell(64,1);
M_.endo_names_tex = cell(64,1);
M_.endo_names_long = cell(64,1);
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
M_.endo_names(13) = {'mps'};
M_.endo_names_tex(13) = {'mps'};
M_.endo_names_long(13) = {'mps'};
M_.endo_names(14) = {'z'};
M_.endo_names_tex(14) = {'z'};
M_.endo_names_long(14) = {'z'};
M_.endo_names(15) = {'EE'};
M_.endo_names_tex(15) = {'EE'};
M_.endo_names_long(15) = {'EE'};
M_.endo_names(16) = {'UE'};
M_.endo_names_tex(16) = {'UE'};
M_.endo_names_long(16) = {'UE'};
M_.endo_names(17) = {'AC'};
M_.endo_names_tex(17) = {'AC'};
M_.endo_names_long(17) = {'AC'};
M_.endo_names(18) = {'C'};
M_.endo_names_tex(18) = {'C'};
M_.endo_names_long(18) = {'C'};
M_.endo_names(19) = {'v'};
M_.endo_names_tex(19) = {'v'};
M_.endo_names_long(19) = {'v'};
M_.endo_names(20) = {'varphi'};
M_.endo_names_tex(20) = {'varphi'};
M_.endo_names_long(20) = {'varphi'};
M_.endo_names(21) = {'VU'};
M_.endo_names_tex(21) = {'VU'};
M_.endo_names_long(21) = {'VU'};
M_.endo_names(22) = {'Q_hat'};
M_.endo_names_tex(22) = {'Q\_hat'};
M_.endo_names_long(22) = {'Q_hat'};
M_.endo_names(23) = {'lambda_hat'};
M_.endo_names_tex(23) = {'lambda\_hat'};
M_.endo_names_long(23) = {'lambda_hat'};
M_.endo_names(24) = {'R_hat'};
M_.endo_names_tex(24) = {'R\_hat'};
M_.endo_names_long(24) = {'R_hat'};
M_.endo_names(25) = {'pi'};
M_.endo_names_tex(25) = {'pi'};
M_.endo_names_long(25) = {'pi'};
M_.endo_names(26) = {'W_hat'};
M_.endo_names_tex(26) = {'W\_hat'};
M_.endo_names_long(26) = {'W_hat'};
M_.endo_names(27) = {'X_hat'};
M_.endo_names_tex(27) = {'X\_hat'};
M_.endo_names_long(27) = {'X_hat'};
M_.endo_names(28) = {'r_u_hat'};
M_.endo_names_tex(28) = {'r\_u\_hat'};
M_.endo_names_long(28) = {'r_u_hat'};
M_.endo_names(29) = {'r_e_hat'};
M_.endo_names_tex(29) = {'r\_e\_hat'};
M_.endo_names_long(29) = {'r_e_hat'};
M_.endo_names(30) = {'mu_e_hat'};
M_.endo_names_tex(30) = {'mu\_e\_hat'};
M_.endo_names_long(30) = {'mu_e_hat'};
M_.endo_names(31) = {'teta_hat'};
M_.endo_names_tex(31) = {'teta\_hat'};
M_.endo_names_long(31) = {'teta_hat'};
M_.endo_names(32) = {'phi_cap_hat'};
M_.endo_names_tex(32) = {'phi\_cap\_hat'};
M_.endo_names_long(32) = {'phi_cap_hat'};
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
M_.endo_names(37) = {'lab_hat'};
M_.endo_names_tex(37) = {'lab\_hat'};
M_.endo_names_long(37) = {'lab_hat'};
M_.endo_names(38) = {'mu_hat'};
M_.endo_names_tex(38) = {'mu\_hat'};
M_.endo_names_long(38) = {'mu_hat'};
M_.endo_names(39) = {'z_hat'};
M_.endo_names_tex(39) = {'z\_hat'};
M_.endo_names_long(39) = {'z_hat'};
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
M_.endo_names(50) = {'g_share'};
M_.endo_names_tex(50) = {'g\_share'};
M_.endo_names_long(50) = {'g_share'};
M_.endo_names(51) = {'H_hat'};
M_.endo_names_tex(51) = {'H\_hat'};
M_.endo_names_long(51) = {'H_hat'};
M_.endo_names(52) = {'TFP_hat'};
M_.endo_names_tex(52) = {'TFP\_hat'};
M_.endo_names_long(52) = {'TFP_hat'};
M_.endo_names(53) = {'EEf_hat'};
M_.endo_names_tex(53) = {'EEf\_hat'};
M_.endo_names_long(53) = {'EEf_hat'};
M_.endo_names(54) = {'UEf_hat'};
M_.endo_names_tex(54) = {'UEf\_hat'};
M_.endo_names_long(54) = {'UEf_hat'};
M_.endo_names(55) = {'AUX_ENDO_LAG_24_1'};
M_.endo_names_tex(55) = {'AUX\_ENDO\_LAG\_24\_1'};
M_.endo_names_long(55) = {'AUX_ENDO_LAG_24_1'};
M_.endo_names(56) = {'AUX_ENDO_LAG_24_2'};
M_.endo_names_tex(56) = {'AUX\_ENDO\_LAG\_24\_2'};
M_.endo_names_long(56) = {'AUX_ENDO_LAG_24_2'};
M_.endo_names(57) = {'AUX_ENDO_LAG_24_3'};
M_.endo_names_tex(57) = {'AUX\_ENDO\_LAG\_24\_3'};
M_.endo_names_long(57) = {'AUX_ENDO_LAG_24_3'};
M_.endo_names(58) = {'AUX_ENDO_LAG_24_4'};
M_.endo_names_tex(58) = {'AUX\_ENDO\_LAG\_24\_4'};
M_.endo_names_long(58) = {'AUX_ENDO_LAG_24_4'};
M_.endo_names(59) = {'AUX_ENDO_LAG_24_5'};
M_.endo_names_tex(59) = {'AUX\_ENDO\_LAG\_24\_5'};
M_.endo_names_long(59) = {'AUX_ENDO_LAG_24_5'};
M_.endo_names(60) = {'AUX_ENDO_LAG_24_6'};
M_.endo_names_tex(60) = {'AUX\_ENDO\_LAG\_24\_6'};
M_.endo_names_long(60) = {'AUX_ENDO_LAG_24_6'};
M_.endo_names(61) = {'AUX_ENDO_LAG_24_7'};
M_.endo_names_tex(61) = {'AUX\_ENDO\_LAG\_24\_7'};
M_.endo_names_long(61) = {'AUX_ENDO_LAG_24_7'};
M_.endo_names(62) = {'AUX_ENDO_LAG_24_8'};
M_.endo_names_tex(62) = {'AUX\_ENDO\_LAG\_24\_8'};
M_.endo_names_long(62) = {'AUX_ENDO_LAG_24_8'};
M_.endo_names(63) = {'AUX_ENDO_LAG_24_9'};
M_.endo_names_tex(63) = {'AUX\_ENDO\_LAG\_24\_9'};
M_.endo_names_long(63) = {'AUX_ENDO_LAG_24_9'};
M_.endo_names(64) = {'AUX_ENDO_LAG_24_10'};
M_.endo_names_tex(64) = {'AUX\_ENDO\_LAG\_24\_10'};
M_.endo_names_long(64) = {'AUX_ENDO_LAG_24_10'};
M_.endo_partitions = struct();
M_.param_names = cell(61,1);
M_.param_names_tex = cell(61,1);
M_.param_names_long = cell(61,1);
M_.param_names(1) = {'s'};
M_.param_names_tex(1) = {'s'};
M_.param_names_long(1) = {'s'};
M_.param_names(2) = {'y_b'};
M_.param_names_tex(2) = {'y\_b'};
M_.param_names_long(2) = {'y_b'};
M_.param_names(3) = {'y_g'};
M_.param_names_tex(3) = {'y\_g'};
M_.param_names_long(3) = {'y_g'};
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
M_.param_names(19) = {'b'};
M_.param_names_tex(19) = {'b'};
M_.param_names_long(19) = {'b'};
M_.param_names(20) = {'eta'};
M_.param_names_tex(20) = {'eta'};
M_.param_names_long(20) = {'eta'};
M_.param_names(21) = {'sz'};
M_.param_names_tex(21) = {'sz'};
M_.param_names_long(21) = {'sz'};
M_.param_names(22) = {'kap_s'};
M_.param_names_tex(22) = {'kap\_s'};
M_.param_names_long(22) = {'kap_s'};
M_.param_names(23) = {'c_p'};
M_.param_names_tex(23) = {'c\_p'};
M_.param_names_long(23) = {'c_p'};
M_.param_names(24) = {'phi_cap_ss'};
M_.param_names_tex(24) = {'phi\_cap\_ss'};
M_.param_names_long(24) = {'phi_cap_ss'};
M_.param_names(25) = {'rho_r'};
M_.param_names_tex(25) = {'rho\_r'};
M_.param_names_long(25) = {'rho_r'};
M_.param_names(26) = {'phi_pi'};
M_.param_names_tex(26) = {'phi\_pi'};
M_.param_names_long(26) = {'phi_pi'};
M_.param_names(27) = {'phi_q'};
M_.param_names_tex(27) = {'phi\_q'};
M_.param_names_long(27) = {'phi_q'};
M_.param_names(28) = {'phi_alt'};
M_.param_names_tex(28) = {'phi\_alt'};
M_.param_names_long(28) = {'phi_alt'};
M_.param_names(29) = {'mu_u'};
M_.param_names_tex(29) = {'mu\_u'};
M_.param_names_long(29) = {'mu_u'};
M_.param_names(30) = {'sigma_z'};
M_.param_names_tex(30) = {'sigma\_z'};
M_.param_names_long(30) = {'sigma_z'};
M_.param_names(31) = {'sigma_mps'};
M_.param_names_tex(31) = {'sigma\_mps'};
M_.param_names_long(31) = {'sigma_mps'};
M_.param_names(32) = {'sigma_mu'};
M_.param_names_tex(32) = {'sigma\_mu'};
M_.param_names_long(32) = {'sigma_mu'};
M_.param_names(33) = {'sigma_phi0'};
M_.param_names_tex(33) = {'sigma\_phi0'};
M_.param_names_long(33) = {'sigma_phi0'};
M_.param_names(34) = {'b_cap'};
M_.param_names_tex(34) = {'b\_cap'};
M_.param_names_long(34) = {'b_cap'};
M_.param_names(35) = {'varb_cap'};
M_.param_names_tex(35) = {'varb\_cap'};
M_.param_names_long(35) = {'varb_cap'};
M_.param_names(36) = {'prox_u'};
M_.param_names_tex(36) = {'prox\_u'};
M_.param_names_long(36) = {'prox_u'};
M_.param_names(37) = {'prox_e'};
M_.param_names_tex(37) = {'prox\_e'};
M_.param_names_long(37) = {'prox_e'};
M_.param_names(38) = {'Q_ss'};
M_.param_names_tex(38) = {'Q\_ss'};
M_.param_names_long(38) = {'Q_ss'};
M_.param_names(39) = {'lambda_ss'};
M_.param_names_tex(39) = {'lambda\_ss'};
M_.param_names_long(39) = {'lambda_ss'};
M_.param_names(40) = {'R_ss'};
M_.param_names_tex(40) = {'R\_ss'};
M_.param_names_long(40) = {'R_ss'};
M_.param_names(41) = {'X_ss'};
M_.param_names_tex(41) = {'X\_ss'};
M_.param_names_long(41) = {'X_ss'};
M_.param_names(42) = {'r_u_ss'};
M_.param_names_tex(42) = {'r\_u\_ss'};
M_.param_names_long(42) = {'r_u_ss'};
M_.param_names(43) = {'r_e_ss'};
M_.param_names_tex(43) = {'r\_e\_ss'};
M_.param_names_long(43) = {'r_e_ss'};
M_.param_names(44) = {'mu_e_ss'};
M_.param_names_tex(44) = {'mu\_e\_ss'};
M_.param_names_long(44) = {'mu_e_ss'};
M_.param_names(45) = {'teta_ss'};
M_.param_names_tex(45) = {'teta\_ss'};
M_.param_names_long(45) = {'teta_ss'};
M_.param_names(46) = {'U_ss'};
M_.param_names_tex(46) = {'U\_ss'};
M_.param_names_long(46) = {'U_ss'};
M_.param_names(47) = {'l_b_ss'};
M_.param_names_tex(47) = {'l\_b\_ss'};
M_.param_names_long(47) = {'l_b_ss'};
M_.param_names(48) = {'l_g_ss'};
M_.param_names_tex(48) = {'l\_g\_ss'};
M_.param_names_long(48) = {'l_g_ss'};
M_.param_names(49) = {'EE_ss'};
M_.param_names_tex(49) = {'EE\_ss'};
M_.param_names_long(49) = {'EE_ss'};
M_.param_names(50) = {'UE_ss'};
M_.param_names_tex(50) = {'UE\_ss'};
M_.param_names_long(50) = {'UE_ss'};
M_.param_names(51) = {'EU_ss'};
M_.param_names_tex(51) = {'EU\_ss'};
M_.param_names_long(51) = {'EU_ss'};
M_.param_names(52) = {'AC_ss'};
M_.param_names_tex(52) = {'AC\_ss'};
M_.param_names_long(52) = {'AC_ss'};
M_.param_names(53) = {'EEf_ss'};
M_.param_names_tex(53) = {'EEf\_ss'};
M_.param_names_long(53) = {'EEf_ss'};
M_.param_names(54) = {'b_share_ss'};
M_.param_names_tex(54) = {'b\_share\_ss'};
M_.param_names_long(54) = {'b_share_ss'};
M_.param_names(55) = {'g_share_ss'};
M_.param_names_tex(55) = {'g\_share\_ss'};
M_.param_names_long(55) = {'g_share_ss'};
M_.param_names(56) = {'v_ss'};
M_.param_names_tex(56) = {'v\_ss'};
M_.param_names_long(56) = {'v_ss'};
M_.param_names(57) = {'C_ss'};
M_.param_names_tex(57) = {'C\_ss'};
M_.param_names_long(57) = {'C_ss'};
M_.param_names(58) = {'varphi_ss'};
M_.param_names_tex(58) = {'varphi\_ss'};
M_.param_names_long(58) = {'varphi_ss'};
M_.param_names(59) = {'empl_ss'};
M_.param_names_tex(59) = {'empl\_ss'};
M_.param_names_long(59) = {'empl_ss'};
M_.param_names(60) = {'H_ss'};
M_.param_names_tex(60) = {'H\_ss'};
M_.param_names_long(60) = {'H_ss'};
M_.param_names(61) = {'TFP_ss'};
M_.param_names_tex(61) = {'TFP\_ss'};
M_.param_names_long(61) = {'TFP_ss'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 64;
M_.param_nbr = 61;
M_.orig_endo_nbr = 54;
M_.aux_vars(1).endo_index = 55;
M_.aux_vars(1).type = 1;
M_.aux_vars(1).orig_index = 25;
M_.aux_vars(1).orig_lead_lag = -1;
M_.aux_vars(1).orig_expr = 'pi(-1)';
M_.aux_vars(2).endo_index = 56;
M_.aux_vars(2).type = 1;
M_.aux_vars(2).orig_index = 25;
M_.aux_vars(2).orig_lead_lag = -2;
M_.aux_vars(2).orig_expr = 'AUX_ENDO_LAG_24_1(-1)';
M_.aux_vars(3).endo_index = 57;
M_.aux_vars(3).type = 1;
M_.aux_vars(3).orig_index = 25;
M_.aux_vars(3).orig_lead_lag = -3;
M_.aux_vars(3).orig_expr = 'AUX_ENDO_LAG_24_2(-1)';
M_.aux_vars(4).endo_index = 58;
M_.aux_vars(4).type = 1;
M_.aux_vars(4).orig_index = 25;
M_.aux_vars(4).orig_lead_lag = -4;
M_.aux_vars(4).orig_expr = 'AUX_ENDO_LAG_24_3(-1)';
M_.aux_vars(5).endo_index = 59;
M_.aux_vars(5).type = 1;
M_.aux_vars(5).orig_index = 25;
M_.aux_vars(5).orig_lead_lag = -5;
M_.aux_vars(5).orig_expr = 'AUX_ENDO_LAG_24_4(-1)';
M_.aux_vars(6).endo_index = 60;
M_.aux_vars(6).type = 1;
M_.aux_vars(6).orig_index = 25;
M_.aux_vars(6).orig_lead_lag = -6;
M_.aux_vars(6).orig_expr = 'AUX_ENDO_LAG_24_5(-1)';
M_.aux_vars(7).endo_index = 61;
M_.aux_vars(7).type = 1;
M_.aux_vars(7).orig_index = 25;
M_.aux_vars(7).orig_lead_lag = -7;
M_.aux_vars(7).orig_expr = 'AUX_ENDO_LAG_24_6(-1)';
M_.aux_vars(8).endo_index = 62;
M_.aux_vars(8).type = 1;
M_.aux_vars(8).orig_index = 25;
M_.aux_vars(8).orig_lead_lag = -8;
M_.aux_vars(8).orig_expr = 'AUX_ENDO_LAG_24_7(-1)';
M_.aux_vars(9).endo_index = 63;
M_.aux_vars(9).type = 1;
M_.aux_vars(9).orig_index = 25;
M_.aux_vars(9).orig_lead_lag = -9;
M_.aux_vars(9).orig_expr = 'AUX_ENDO_LAG_24_8(-1)';
M_.aux_vars(10).endo_index = 64;
M_.aux_vars(10).type = 1;
M_.aux_vars(10).orig_index = 25;
M_.aux_vars(10).orig_lead_lag = -10;
M_.aux_vars(10).orig_expr = 'AUX_ENDO_LAG_24_9(-1)';
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
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
M_.orig_eq_nbr = 54;
M_.eq_nbr = 64;
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
 0 19 0;
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
 0 41 83;
 1 42 0;
 2 43 84;
 0 44 85;
 0 45 0;
 0 46 0;
 0 47 0;
 3 48 0;
 0 49 0;
 0 50 0;
 4 51 0;
 5 52 0;
 0 53 0;
 0 54 0;
 0 55 0;
 6 56 0;
 7 57 0;
 8 58 0;
 0 59 0;
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
 9 73 0;
 10 74 0;
 11 75 0;
 12 76 0;
 13 77 0;
 14 78 0;
 15 79 0;
 16 80 0;
 17 81 0;
 18 82 0;]';
M_.nstatic = 44;
M_.nfwrd   = 2;
M_.npred   = 17;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 18;
M_.ndynamic   = 20;
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
  12 , 'name' , 'phi_cap_hat' ;
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
  28 , 'name' , 'g_share' ;
  29 , 'name' , 'H_hat' ;
  30 , 'name' , 'lab_hat' ;
  31 , 'name' , 'TFP_hat' ;
  32 , 'name' , 'EEf_hat' ;
  33 , 'name' , 'UEf_hat' ;
  34 , 'name' , 'Q' ;
  35 , 'name' , 'lambda' ;
  36 , 'name' , 'R' ;
  37 , 'name' , 'X' ;
  38 , 'name' , 'z' ;
  39 , 'name' , 'r_u' ;
  40 , 'name' , 'r_e' ;
  41 , 'name' , 'mu' ;
  42 , 'name' , 'teta' ;
  43 , 'name' , 'U' ;
  44 , 'name' , 'l_b' ;
  45 , 'name' , 'l_g' ;
  46 , 'name' , 'mps' ;
  47 , 'name' , 'EE' ;
  48 , 'name' , 'UE' ;
  49 , 'name' , 'AC' ;
  50 , 'name' , 'C' ;
  51 , 'name' , 'v' ;
  52 , 'name' , 'varphi' ;
  53 , 'name' , 'empl' ;
  54 , 'name' , 'VU' ;
};
M_.mapping.Q.eqidx = [34 ];
M_.mapping.lambda.eqidx = [35 54 ];
M_.mapping.R.eqidx = [36 ];
M_.mapping.X.eqidx = [37 ];
M_.mapping.r_u.eqidx = [39 ];
M_.mapping.r_e.eqidx = [40 ];
M_.mapping.teta.eqidx = [42 ];
M_.mapping.U.eqidx = [43 ];
M_.mapping.l_b.eqidx = [44 ];
M_.mapping.l_g.eqidx = [45 ];
M_.mapping.empl.eqidx = [53 ];
M_.mapping.mu.eqidx = [41 ];
M_.mapping.mps.eqidx = [46 ];
M_.mapping.z.eqidx = [38 ];
M_.mapping.EE.eqidx = [47 ];
M_.mapping.UE.eqidx = [48 ];
M_.mapping.AC.eqidx = [49 ];
M_.mapping.C.eqidx = [50 ];
M_.mapping.v.eqidx = [51 ];
M_.mapping.varphi.eqidx = [52 ];
M_.mapping.VU.eqidx = [54 ];
M_.mapping.Q_hat.eqidx = [6 11 14 15 31 34 ];
M_.mapping.lambda_hat.eqidx = [1 2 3 4 6 14 23 26 29 35 ];
M_.mapping.R_hat.eqidx = [1 15 36 ];
M_.mapping.pi.eqidx = [1 11 15 19 ];
M_.mapping.W_hat.eqidx = [2 3 4 26 ];
M_.mapping.X_hat.eqidx = [2 3 4 11 13 14 25 29 37 ];
M_.mapping.r_u_hat.eqidx = [3 7 9 13 20 21 39 ];
M_.mapping.r_e_hat.eqidx = [4 7 13 20 40 ];
M_.mapping.mu_e_hat.eqidx = [4 5 ];
M_.mapping.teta_hat.eqidx = [6 7 9 12 13 20 21 24 42 ];
M_.mapping.phi_cap_hat.eqidx = [12 ];
M_.mapping.U_hat.eqidx = [6 7 8 9 10 13 14 20 24 33 43 ];
M_.mapping.l_b_hat.eqidx = [5 7 8 14 20 27 44 ];
M_.mapping.l_g_hat.eqidx = [8 28 45 ];
M_.mapping.empl_hat.eqidx = [5 10 27 28 30 31 32 53 ];
M_.mapping.lab_hat.eqidx = [30 ];
M_.mapping.mu_hat.eqidx = [6 14 17 23 41 ];
M_.mapping.z_hat.eqidx = [2 3 4 13 14 16 25 29 38 ];
M_.mapping.mps_hat.eqidx = [15 18 46 ];
M_.mapping.INFL.eqidx = [19 ];
M_.mapping.EE_hat.eqidx = [20 22 32 47 ];
M_.mapping.UE_hat.eqidx = [21 22 33 48 ];
M_.mapping.AC_hat.eqidx = [15 22 49 ];
M_.mapping.C_hat.eqidx = [23 50 ];
M_.mapping.v_hat.eqidx = [24 51 ];
M_.mapping.varphi_hat.eqidx = [25 52 ];
M_.mapping.WplusLam.eqidx = [26 ];
M_.mapping.b_share.eqidx = [27 ];
M_.mapping.g_share.eqidx = [28 ];
M_.mapping.H_hat.eqidx = [29 30 31 ];
M_.mapping.TFP_hat.eqidx = [31 ];
M_.mapping.EEf_hat.eqidx = [32 ];
M_.mapping.UEf_hat.eqidx = [33 ];
M_.mapping.eps_mu.eqidx = [17 ];
M_.mapping.eps_z.eqidx = [16 ];
M_.mapping.eps_mps.eqidx = [18 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [24 25 30 33 34 38 39 40 55 56 57 58 59 60 61 62 63 64 ];
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(64, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(61, 1);
M_.endo_trends = struct('deflator', cell(64, 1), 'log_deflator', cell(64, 1), 'growth_factor', cell(64, 1), 'log_growth_factor', cell(64, 1));
M_.NNZDerivatives = [202; -1; -1; ];
M_.static_tmp_nbr = [15; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 0.5;
s = M_.params(1);
M_.params(2) = 1;
y_b = M_.params(2);
M_.params(3) = 1.8;
y_g = M_.params(3);
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
M_.params(49) = 0.024;
EE_ss = M_.params(49);
M_.params(50) = 0.25;
UE_ss = M_.params(50);
M_.params(51) = 0.01314;
EU_ss = M_.params(51);
M_.params(52) = M_.params(49)/M_.params(50);
AC_ss = M_.params(52);
M_.params(46) = M_.params(51)/(M_.params(50)+M_.params(51));
U_ss = M_.params(46);
M_.params(59) = 1-M_.params(46);
empl_ss = M_.params(59);
M_.params(53) = M_.params(49)*M_.params(59);
EEf_ss = M_.params(53);
M_.params(45) = 1;
teta_ss = M_.params(45);
M_.params(24) = 1;
phi_cap_ss = M_.params(24);
M_.params(20) = 6;
eta = M_.params(20);
M_.params(17) = 0.995;
beta = M_.params(17);
M_.params(12) = 0.94;
rho_mu = M_.params(12);
M_.params(14) = 0.494;
rho_z = M_.params(14);
M_.params(15) = 0.5;
rho_phi0 = M_.params(15);
M_.params(40) = 1/M_.params(17)-1;
R_ss = M_.params(40);
M_.params(30) = 0.3;
sigma_z = M_.params(30);
M_.params(31) = 0.00005;
sigma_mps = M_.params(31);
M_.params(32) = 0.01;
sigma_mu = M_.params(32);
M_.params(13) = 0.93;
rho_mps = M_.params(13);
M_.params(25) = 0.85;
rho_r = M_.params(25);
M_.params(18) = 1-M_.params(4);
xi_b = M_.params(18);
M_.params(16) = M_.params(4)*(M_.params(49)+M_.params(51))/(M_.params(4)+M_.params(4)*M_.params(18)-M_.params(52));
delta = M_.params(16);
M_.params(21) = (1-M_.params(51)/M_.params(16))/M_.params(50);
sz = M_.params(21);
r_ratio = (s*AC_ss*(1-delta)*UE_ss)/(EE_ss-delta+EU_ss); 
M_.params(47) = M_.params(18)*(1-M_.params(46))*M_.params(16)/(M_.params(16)+M_.params(4)*(M_.params(49)-M_.params(50)*M_.params(16)*M_.params(21))/M_.params(52));
l_b_ss = M_.params(47);
M_.params(48) = 1-M_.params(46)-M_.params(47);
l_g_ss = M_.params(48);
M_.params(54) = M_.params(47)/M_.params(59);
b_share_ss = M_.params(54);
M_.params(55) = M_.params(48)/M_.params(59);
g_share_ss = M_.params(55);
M_.params(29) = M_.params(4)*M_.params(3)^(1+M_.params(10))+M_.params(18)*M_.params(2)^(1+M_.params(10));
mu_u = M_.params(29);
M_.params(44) = M_.params(4)*M_.params(47)/(1-M_.params(46))*(M_.params(3)^(1+M_.params(10))-M_.params(2)^(1+M_.params(10)));
mu_e_ss = M_.params(44);
M_.params(39) = lambda_solver(M_.params(5),M_.params(6),M_.params(20),M_.params(17),M_.params(7),M_.params(10),M_.params(8),M_.params(1),M_.params(47),M_.params(48),M_.params(44),M_.params(2),M_.params(3),M_.params(16),M_.params(21),M_.params(49),M_.params(50),M_.params(51),M_.params(52),M_.params(46));
lambda_ss = M_.params(39);
M_.params(38) = Q_solver(M_.params(5),M_.params(6),M_.params(20),M_.params(17),M_.params(7),M_.params(10),M_.params(8),M_.params(1),M_.params(47),M_.params(48),M_.params(44),M_.params(2),M_.params(3),M_.params(16),M_.params(21),M_.params(49),M_.params(50),M_.params(51),M_.params(52),M_.params(46));
Q_ss = M_.params(38);
M_.params(41) = M_.params(5)*(M_.params(20)-1)/M_.params(20)*M_.params(38)^((M_.params(5)-1)/M_.params(5));
X_ss = M_.params(41);
M_.params(35) = 1/(1+M_.params(10))*(M_.params(41)*M_.params(39)/M_.params(8))^M_.params(10);
varb_cap = M_.params(35);
M_.params(34) = M_.params(35)*M_.params(29)-M_.params(44)*M_.params(35)*(M_.params(50)*M_.params(52)*M_.params(1)*(1-M_.params(16))/(M_.params(51)+M_.params(49)-M_.params(16)))^M_.params(7);
b_cap = M_.params(34);
M_.params(19) = M_.params(39)*M_.params(41)*M_.params(34);
b = M_.params(19);
M_.params(36) = M_.params(17)/(1-M_.params(17)*(1-M_.params(16)))*(M_.params(35)*M_.params(29)-M_.params(34));
prox_u = M_.params(36);
M_.params(37) = M_.params(44)*M_.params(35)*M_.params(17)/(1-M_.params(17)*(1-M_.params(16)));
prox_e = M_.params(37);
M_.params(23) = (M_.params(38)-M_.params(39)^(-M_.params(6)))/(M_.params(46)+(1-M_.params(46))*(M_.params(16)*M_.params(21)+M_.params(1)*(1-M_.params(16))));
c_p = M_.params(23);
M_.params(42) = M_.params(50);
r_u_ss = M_.params(42);
M_.params(43) = (M_.params(49)-M_.params(50)*M_.params(16)*M_.params(21))/(M_.params(1)*M_.params(52)*(1-M_.params(16)));
r_e_ss = M_.params(43);
M_.params(22) = M_.params(36)/M_.params(42)^M_.params(7);
kap_s = M_.params(22);
M_.params(57) = M_.params(39)^(-M_.params(6));
C_ss = M_.params(57);
M_.params(56) = M_.params(46)+(1-M_.params(46))*(M_.params(16)*M_.params(21)+M_.params(1)*(1-M_.params(16)));
v_ss = M_.params(56);
M_.params(58) = M_.params(41);
varphi_ss = M_.params(58);
VU_ss = (beta*b)/(lambda_ss*(1-beta*(1-delta)));
M_.params(60) = (M_.params(41)*M_.params(39)/M_.params(8))^M_.params(10)*(M_.params(47)*M_.params(2)^M_.params(10)+(M_.params(48)-M_.params(47))*M_.params(3)^M_.params(10))/(1-M_.params(46));
H_ss = M_.params(60);
M_.params(61) = M_.params(38)/M_.params(59)*M_.params(60);
TFP_ss = M_.params(61);
M_.params(3) = 1.01;
y_g = M_.params(3);
M_.params(1) = 0.5;
s = M_.params(1);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(32))^2;
M_.Sigma_e(2, 2) = (M_.params(30))^2;
M_.Sigma_e(3, 3) = (M_.params(31))^2;
%
% OPTIM_WEIGHTS
%
M_.osr.variable_weights = sparse(M_.endo_nbr,M_.endo_nbr);
M_.osr.variable_indices = [];

M_.osr.variable_weights(22,22) = 0.25;
M_.osr.variable_indices = [M_.osr.variable_indices; 22];
M_.osr.variable_weights(25,25) = 1;
M_.osr.variable_indices = [M_.osr.variable_indices; 25];
M_.osr.param_names = {'phi_pi';'phi_q';'phi_alt'};
M_.osr.param_names = cellstr(M_.osr.param_names);
M_.osr.param_indices = zeros(length(M_.osr.param_names), 1);
M_.osr.param_indices(1) = 26;
M_.osr.param_indices(2) = 27;
M_.osr.param_indices(3) = 28;
M_.params(26) = 1.17;
phi_pi = M_.params(26);
M_.params(27) = 0;
phi_q = M_.params(27);
M_.params(28) = 0;
phi_alt = M_.params(28);
options_.irf = 0;
options_.osr.opt_algo = 2;
var_list_ = {};
oo_.osr = osr(var_list_,M_.osr.param_names,M_.osr.variable_indices,M_.osr.variable_weights);
steady;
resid;
oo_.dr.eigval = check(M_,options_,oo_);
options_.nograph = 1;   
options_.impulse_responses.plot_threshold = 0;
options_.irf = 30;
options_.order = 1;
options_.periods = 10000;
var_list_ = {'R_hat';'Q_hat';'C_hat';'pi';'lab_hat';'W_hat';'TFP_hat';'U_hat';'H_hat';'b_share';'g_share';'EEf_hat';'UEf_hat'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
options_.impulse_responses.plot_threshold = 0;
options_.irf = 30;
options_.order = 1;
options_.periods = 0;
options_.conditional_variance_decomposition = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;];
var_list_ = {'R_hat';'Q_hat';'C_hat';'pi';'lab_hat';'W_hat';'TFP_hat';'U_hat';'H_hat';'b_share';'g_share';'EEf_hat';'UEf_hat'};
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
disp('Note: 33 warning(s) encountered in the preprocessor')
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
