 %% Housekeeping

clc
clear
set(0,'DefaultFigureWindowStyle','docked'); % I prefer docked figures...

%% Run four different simulations
%% Simulate the benchmark economy

dynare linear_sims -DBenchmarkModel

SavedSims.BenchmarkModel = oo_.endo_simul;
SavedInfo.endo_names.BenchmarkModel  = M_.endo_names;
save SimulationFile SavedSims SavedInfo

%% Simulate the Ramsey optimal economy

dynare linear_sims -DRamseyModel

load SimulationFile
SavedSims.RamseyModel = oo_.endo_simul;
save SimulationFile SavedSims

%% Simulate EE simple rule
dynare linear_sims -DEE_SimpleRule

load SimulationFile
SavedSims.EE_SimpleRule = oo_.endo_simul;
save SimulationFile SavedSims

%% Simulate b_share simple rule
dynare linear_sims -Db_shareRule

load SimulationFile
SavedSims.b_shareSimpleRule = oo_.endo_simul;
save SimulationFile SavedSims

%% Saving TS for variables in loss function
Y_benchmark = SavedSims.BenchmarkModel(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'Q_hat')),1000:10000);
PI_benchmark = SavedSims.BenchmarkModel(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'INFL')),1000:10000);
U_benchmark = SavedSims.BenchmarkModel(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'U_hat')),1000:10000);

Y_ramsey = SavedSims.RamseyModel(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'Q_hat')),1000:10000);
PI_ramsey = SavedSims.RamseyModel(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'INFL')),1000:10000);
U_ramsey = SavedSims.RamseyModel(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'U_hat')),1000:10000);

Y_EEsr = SavedSims.EE_SimpleRule(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'Q_hat')),1000:10000);
PI_EEsr = SavedSims.EE_SimpleRule(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'INFL')),1000:10000);
U_EEsr = SavedSims.EE_SimpleRule(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'U_hat')),1000:10000);

Y_b_sharesr = SavedSims.b_shareSimpleRule(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'Q_hat')),1000:10000);
PI_b_sharesr = SavedSims.b_shareSimpleRule(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'INFL')),1000:10000);
U_b_sharesr = SavedSims.b_shareSimpleRule(find(strcmp(SavedInfo.endo_names.BenchmarkModel, 'U_hat')),1000:10000);

%% Evaluate Loss
%loss function with output: F = 0.25*(Y^2) + pi^2
EvaluatedLoss.Output.BenchmarkModel = sum(Y_benchmark.^2 * 0.25+ PI_benchmark.^2);
EvaluatedLoss.Output.RamseyModel = sum(Y_ramsey.^2 * 0.25+ PI_ramsey.^2);
EvaluatedLoss.Output.EE_SimpleRule = sum(Y_EEsr.^2 * 0.25+ PI_EEsr.^2);
EvaluatedLoss.Output.b_shareSimpleRule = sum(Y_b_sharesr.^2 * 0.25+ PI_b_sharesr.^2);

%loss function with unemployment: F = (U^2) + pi^2
EvaluatedLoss.Unemployment.BenchmarkModel = sum(U_benchmark.^2 + PI_benchmark.^2);
EvaluatedLoss.Unemployment.RamseyModel = sum(U_ramsey.^2 + PI_ramsey.^2);
EvaluatedLoss.Unemployment.EE_SimpleRule = sum(U_EEsr.^2 + PI_EEsr.^2);
EvaluatedLoss.Unemployment.b_shareSimpleRule= sum(U_b_sharesr.^2 + PI_b_sharesr.^2);
%% Comparison
%e.g. A result of -99 would indicate that by using the specified rule instead
%of the ramsey optimal we are losing 99% of "welfare" created by Ramsey 
% OR Ramsey improves welfare by 99%. 

Ramsey_vs_benchmark = (EvaluatedLoss.Output.RamseyModel/EvaluatedLoss.Output.BenchmarkModel-1)*100;
Ramsey_vs_EEsr = (EvaluatedLoss.Output.RamseyModel/EvaluatedLoss.Output.EE_SimpleRule-1)*100;
Ramsey_vs_b_sharesr = (EvaluatedLoss.Output.RamseyModel/EvaluatedLoss.Output.b_shareSimpleRule-1)*100;

EEsr_vs_benchmark = (EvaluatedLoss.Output.EE_SimpleRule/EvaluatedLoss.Output.BenchmarkModel-1)*100;
b_sharesr_vs_benchmark = (EvaluatedLoss.Output.b_shareSimpleRule/EvaluatedLoss.Output.BenchmarkModel-1)*100;
b_sharesr_vs_EEsr = (EvaluatedLoss.Output.b_shareSimpleRule/EvaluatedLoss.Output.EE_SimpleRule-1)*100;

