 %% Housekeeping

clc
clear
set(0,'DefaultFigureWindowStyle','docked'); % I prefer docked figures...
%% Overview of the structure
%{
Loss function using Output (Q)
    Ramsey model
    All Taylor Rules based on Unemployment (including Benchamark)
    All Taylor Rules based on Output
    Comupute losses compared to Ramsey
    Display table of Taylor Rules
    IRFs
    FEVD

Loss function using Unemployment (U)
    Ramsey model
    All Taylor Rules based on Unemployment (including Benchamark)
    All Taylor Rules based on Output
    Comupute losses compared to Ramsey
    Display table of Taylor Rules
    IRFs
    FEVD        
%}
%% Inital settings

Loss_Q = 0; %if 1: Loss function is computed using Output; if 0: Loss function is computed using Unemployment
Q_Rule = 1;
U_Rule = 1;







%% Loss is Output
if Loss_Q==1
% inital settings
wy = 0.25;
%% Simulate the Ramsey optimal economy 
dynare linear_sims -DRamseyModel -DLoss_Q

SavedInfo.endo_names  = M_.endo_names;
SavedSims.RamseyModel = oo_.endo_simul;
SavedIRF.RamseyModel = oo_.irfs;
SavedFEVD.RamseyModel = oo_.variance_decomposition; 
Saved_conditional_FEVD.RamseyModel = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims

%% Run models with Q rule
if U_Rule==1

%% Simulate the Benchmark economy 

dynare linear_sims -DBenchmarkModel -DLoss_Q
%add SavedInfo
SavedSims.BenchmarkModel = oo_.endo_simul;
SavedOSR.BenchmarkModel.phi_pi = M_.params(find(strcmp(M_.param_names, 'phi_pi')),1);
SavedOSR.BenchmarkModel.phi_u = M_.params(find(strcmp(M_.param_names, 'phi_u')),1);
SavedIRF.BenchmarkModel = oo_.irfs;
SavedFEVD.BenchmarkModel = oo_.variance_decomposition; 
Saved_conditional_FEVD.BenchmarkModel = oo_.conditional_variance_decomposition;
save SimulationFile SavedSims SavedInfo

%% Simulate the benchmark economy with U in TR - Optimized

dynare linear_sims -DU_OptimalTR -DLoss_Q

load SimulationFile
SavedSims.U_OptimalTR = oo_.endo_simul;
SavedOSR.U_OptimalTR = oo_.osr.optim_params;
SavedIRF.U_OptimalTR = oo_.irfs;
SavedFEVD.U_OptimalTR = oo_.variance_decomposition; 
Saved_conditional_FEVD.U_OptimalTR = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims SavedInfo

%% Simulate EE simple rule
dynare linear_sims -DU_EEosr -DLoss_Q


load SimulationFile
SavedSims.U_EEosr = oo_.endo_simul;
SavedOSR.U_EEosr = oo_.osr.optim_params;
SavedIRF.U_EEosr = oo_.irfs;
SavedFEVD.U_EEosr = oo_.variance_decomposition; 
Saved_conditional_FEVD.U_EEosr = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims

%% Simulate Alternative simple rule
dynare linear_sims -DU_ALTosr -DLoss_Q


load SimulationFile
SavedSims.U_ALTosr = oo_.endo_simul;
SavedOSR.U_ALTosr = oo_.osr.optim_params;
SavedIRF.U_ALTosr = oo_.irfs;
SavedFEVD.U_ALTosr = oo_.variance_decomposition; 
Saved_conditional_FEVD.U_ALTosr = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims

end

if Q_Rule==1

%% Simulate the Optimized Taylor Rule - Output

dynare linear_sims -DQ_OptimalTR -DLoss_Q


SavedSims.Q_OptimalTR = oo_.endo_simul;
SavedOSR.Q_OptimalTR = oo_.osr.optim_params;
SavedIRF.Q_OptimalTR = oo_.irfs;
SavedFEVD.Q_OptimalTR = oo_.variance_decomposition; 
Saved_conditional_FEVD.Q_OptimalTR = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims SavedInfo


%% Simulate EE simple rule
dynare linear_sims -DQ_EEosr -DLoss_Q


load SimulationFile
SavedSims.Q_EEosr = oo_.endo_simul;
SavedOSR.Q_EEosr = oo_.osr.optim_params;
SavedIRF.Q_EEosr = oo_.irfs;
SavedFEVD.Q_EEosr = oo_.variance_decomposition; 
Saved_conditional_FEVD.Q_EEosr = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims

%% Simulate Alternative simple rule
dynare linear_sims -DQ_ALTosr -DLoss_Q


load SimulationFile
SavedSims.Q_ALTosr = oo_.endo_simul;
SavedOSR.Q_ALTosr = oo_.osr.optim_params;
SavedIRF.Q_ALTosr = oo_.irfs;
SavedFEVD.Q_ALTosr = oo_.variance_decomposition; 
Saved_conditional_FEVD.Q_ALTosr = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims

end

%% End in running model
%% Automatic way to save variables for loss function and compute loss function

fields = fieldnames(SavedSims);
fieldsToRemove = {'RamseyModel'}; %removing these fields 
fields = fields(~ismember(fields, fieldsToRemove));

%compte Loss for Ramsey
Q_ramsey = SavedSims.RamseyModel(find(strcmp(SavedInfo.endo_names, 'Q_hat')),1000:10000);
PI_ramsey = SavedSims.RamseyModel(find(strcmp(SavedInfo.endo_names, 'Q_hat')),1000:10000);
LOSS_ramsey = sum(Q_ramsey.^2 * wy+ PI_ramsey.^2);

%compute Loss to Ramsey for other Models
for i=1:length(fields)
    Q = SavedSims.(fields{i})(find(strcmp(SavedInfo.endo_names, 'Q_hat')),1000:10000);
    PI = SavedSims.(fields{i})(find(strcmp(SavedInfo.endo_names, 'pi')),1000:10000);

    LOSS = sum(Q.^2 * wy+ PI.^2);
    Loss_to_Ramsey.Output.(fields{i}) = (LOSS_ramsey/LOSS-1)*100;
end

%% Plot table 
allVariables = ["phi_pi","phi_u","phi_q","phi_EE", "phi_alt"];
numVariables = length(allVariables);
dataTable = array2table(nan(length(fields), numVariables), 'RowNames', fields, 'VariableNames', allVariables);

% Step 3: Populate the table
for i = 1:length(fields)
    for j = 1:numVariables
        variableName = allVariables{j};
        if isfield(SavedOSR.(fields{i}), variableName)
            dataTable{i, variableName} = SavedOSR.(fields{i}).(variableName);
        end
    end
end

lossValues = struct2array(Loss_to_Ramsey.Output);
lossTable = array2table(lossValues, 'VariableNames', {'Loss to Ramsey'});
% Concatenate the new table with the existing dataTable
dataTable = [lossTable, dataTable];

% Display the table
disp(dataTable);

%% IRF production 
% Assuming SavedIRF is your main structure
fields = fieldnames(SavedIRF); % Get the names of the models
irfVars = fieldnames(SavedIRF.(fields{1})); % Get IRF variable names from the first model
outputDir = fullfile(pwd, 'IRFs'); % Output directory

for i = 1:length(irfVars)
    f=figure; % Create a new figure for each IRF variable
    hold on; % Hold on to plot multiple lines in the same figure
    for j = 1:length(fields)
        % Extracting the IRF data for each model
        irfData = SavedIRF.(fields{j}).(irfVars{i});
        legend_name = strrep(fields{j}, '_', ' ');
        plot(irfData, 'DisplayName', legend_name, 'LineWidth', 2); % Plot with legend label
    end
    hold off;
    legend('show'); % Show legend
    formattedVarName = strrep(irfVars{i}, '_', ' ');
    title(['IRF for ' formattedVarName]); % Title for each plot
    xlabel('Periods'); % X-axis label
    ylabel('Response'); % Y-axis label

    saveas(f, fullfile(outputDir, ['IRF_for_' formattedVarName '.png']));
end

%% saving the FEVD - unconditional variance 
var_names = oo_.var_list;
shock_names = (M_.exo_names)';
fields = fieldnames(SavedFEVD); % Get the names of the models

for i = 1:length(fields)
    T_fevd.(fields{i}) = array2table(SavedFEVD.(fields{i}), 'VariableNames', shock_names, 'RowNames', var_names);
end

disp(T_fevd.BenchmarkModel);

end












%% Loss is Unemployment
if Loss_Q==0
% inital settings
wy = 0.25;
%% Simulate the Ramsey optimal economy 
dynare linear_sims -DRamseyModel -DLoss_Q

SavedInfo.endo_names  = M_.endo_names;
SavedSims.RamseyModel = oo_.endo_simul;
SavedIRF.RamseyModel = oo_.irfs;
SavedFEVD.RamseyModel = oo_.variance_decomposition; 
Saved_conditional_FEVD.RamseyModel = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims

%% Run models with Q rule
if U_Rule==1

%% Simulate the Benchmark economy 

dynare linear_sims -DBenchmarkModel -DLoss_Q
%add SavedInfo
SavedSims.BenchmarkModel = oo_.endo_simul;
SavedOSR.BenchmarkModel.phi_pi = M_.params(find(strcmp(M_.param_names, 'phi_pi')),1);
SavedOSR.BenchmarkModel.phi_u = M_.params(find(strcmp(M_.param_names, 'phi_u')),1);
SavedIRF.BenchmarkModel = oo_.irfs;
SavedFEVD.BenchmarkModel = oo_.variance_decomposition; 
Saved_conditional_FEVD.BenchmarkModel = oo_.conditional_variance_decomposition;
save SimulationFile SavedSims SavedInfo

%% Simulate the benchmark economy with U in TR - Optimized

dynare linear_sims -DU_OptimalTR -DLoss_Q


load SimulationFile
SavedSims.U_OptimalTR = oo_.endo_simul;
SavedOSR.U_OptimalTR = oo_.osr.optim_params;
SavedIRF.U_OptimalTR = oo_.irfs;
SavedFEVD.U_OptimalTR = oo_.variance_decomposition; 
Saved_conditional_FEVD.U_OptimalTR = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims SavedInfo

%% Simulate EE simple rule
dynare linear_sims -DU_EEosr -DLoss_Q


load SimulationFile
SavedSims.U_EEosr = oo_.endo_simul;
SavedOSR.U_EEosr = oo_.osr.optim_params;
SavedIRF.U_EEosr = oo_.irfs;
SavedFEVD.U_EEosr = oo_.variance_decomposition; 
Saved_conditional_FEVD.U_EEosr = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims

%% Simulate Alternative simple rule
dynare linear_sims -DU_ALTosr -DLoss_Q


load SimulationFile
SavedSims.U_ALTosr = oo_.endo_simul;
SavedOSR.U_ALTosr = oo_.osr.optim_params;
SavedIRF.U_ALTosr = oo_.irfs;
SavedFEVD.U_ALTosr = oo_.variance_decomposition; 
Saved_conditional_FEVD.U_ALTosr = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims

end

if Q_Rule==1

%% Simulate the Optimized Taylor Rule - Output

dynare linear_sims -DQ_OptimalTR -DLoss_Q


SavedSims.Q_OptimalTR = oo_.endo_simul;
SavedOSR.Q_OptimalTR = oo_.osr.optim_params;
SavedIRF.Q_OptimalTR = oo_.irfs;
SavedFEVD.Q_OptimalTR = oo_.variance_decomposition; 
Saved_conditional_FEVD.Q_OptimalTR = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims SavedInfo


%% Simulate EE simple rule
dynare linear_sims -DQ_EEosr -DLoss_Q


load SimulationFile
SavedSims.Q_EEosr = oo_.endo_simul;
SavedOSR.Q_EEosr = oo_.osr.optim_params;
SavedIRF.Q_EEosr = oo_.irfs;
SavedFEVD.Q_EEosr = oo_.variance_decomposition; 
Saved_conditional_FEVD.Q_EEosr = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims

%% Simulate Alternative simple rule
dynare linear_sims -DQ_ALTosr -DLoss_Q


load SimulationFile
SavedSims.Q_ALTosr = oo_.endo_simul;
SavedOSR.Q_ALTosr = oo_.osr.optim_params;
SavedIRF.Q_ALTosr = oo_.irfs;
SavedFEVD.Q_ALTosr = oo_.variance_decomposition; 
Saved_conditional_FEVD.Q_ALTosr = oo_.conditional_variance_decomposition;

save SimulationFile SavedSims

end

%% End in running model
%% Automatic way to save variables for loss function and compute loss function

fields = fieldnames(SavedSims);
fieldsToRemove = { 'RamseyModel'}; %removing these fields 'BenchmarkModel_no_optim'
fields = fields(~ismember(fields, fieldsToRemove));

%compte Loss for Ramsey
Q_ramsey = SavedSims.RamseyModel(find(strcmp(SavedInfo.endo_names, 'Q_hat')),1000:10000);
PI_ramsey = SavedSims.RamseyModel(find(strcmp(SavedInfo.endo_names, 'Q_hat')),1000:10000);
LOSS_ramsey = sum(Q_ramsey.^2 * wy+ PI_ramsey.^2);

%compute Loss to Ramsey for other Models
for i=1:length(fields)
    Q = SavedSims.(fields{i})(find(strcmp(SavedInfo.endo_names, 'Q_hat')),1000:10000);
    PI = SavedSims.(fields{i})(find(strcmp(SavedInfo.endo_names, 'pi')),1000:10000);

    LOSS = sum(Q.^2 * wy+ PI.^2);
    Loss_to_Ramsey.Unemployment.(fields{i}) = (LOSS_ramsey/LOSS-1)*100;
end

%% Plot table 
allVariables = ["phi_pi","phi_u","phi_q","phi_EE", "phi_alt"];
numVariables = length(allVariables);
dataTable = array2table(nan(length(fields), numVariables), 'RowNames', fields, 'VariableNames', allVariables);

% Step 3: Populate the table
for i = 1:length(fields)
    for j = 1:numVariables
        variableName = allVariables{j};
        if isfield(SavedOSR.(fields{i}), variableName)
            dataTable{i, variableName} = SavedOSR.(fields{i}).(variableName);
        end
    end
end

lossValues = struct2array(Loss_to_Ramsey.Unemployment);
lossTable = array2table(lossValues, 'VariableNames', {'Loss to Ramsey'});
% Concatenate the new table with the existing dataTable
dataTable = [lossTable, dataTable];

% Display the table
disp(dataTable);

%% IRF production 
% Assuming SavedIRF is your main structure
fields = fieldnames(SavedIRF); % Get the names of the models
irfVars = fieldnames(SavedIRF.(fields{1})); % Get IRF variable names from the first model
outputDir = fullfile(pwd, 'IRFs'); % Output directory

for i = 1:length(irfVars)
    f=figure; % Create a new figure for each IRF variable
    hold on; % Hold on to plot multiple lines in the same figure
    for j = 1:length(fields)
        % Extracting the IRF data for each model
        irfData = SavedIRF.(fields{j}).(irfVars{i});
        legend_name = strrep(fields{j}, '_', ' ');
        plot(irfData, 'DisplayName', legend_name, 'LineWidth', 2); % Plot with legend label
    end
    hold off;
    legend('show'); % Show legend
    formattedVarName = strrep(irfVars{i}, '_', ' ');
    title(['IRF for ' formattedVarName]); % Title for each plot
    xlabel('Periods'); % X-axis label
    ylabel('Response'); % Y-axis label

    saveas(f, fullfile(outputDir, ['IRF_for_' formattedVarName '.png']));
end

%% saving the FEVD - unconditional variance 
var_names = oo_.var_list;
shock_names = (M_.exo_names)';
fields = fieldnames(SavedFEVD); % Get the names of the models

for i = 1:length(fields)
    T_fevd.(fields{i}) = array2table(SavedFEVD.(fields{i}), 'VariableNames', shock_names, 'RowNames', var_names);
end

disp(T_fevd.BenchmarkModel);

end




