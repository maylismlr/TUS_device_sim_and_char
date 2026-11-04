%% Plot delay calibration and create validation table
%
% Loads the previously computed delay-calibration data after running get_validation_delays.m.
%
% Outputs:
%   1. Inverse PCHIP calibration figure
%   2. Validation results table
%
% Plots Supp. Figue 6.5 and partial results from Table 3.3 from the
% Master's thesis.
% 
% author: Maylis Mueller
% date:
% last update:

clear;
clc;
close all;

save_outputs = false;

% =========================================================================
% LOAD DATA
% =========================================================================

% Calibration data: contains the PCHIP-derived corrected inputs
cal = load('delay_calib_results/delay_calibration.mat');

% Validation data: contains before/after simulation results
val = load('delay_calib_results/delay_calibration_validation.mat');

% Reference focal settings
focus_lengths = val.focus_lengths;
tus_device    = val.tus_device;

% Uncorrected simulated peak positions
our_sim = tus_device + val.uncorrected_error;

% Corrected simulation inputs from PCHIP calibration
corrected_sim_input = val.corrected_input;

% Corrected simulated peak positions
corrected_sim_peak = val.corrected_sim_peak;

% COLORS
orange = [0.8500 0.3250 0.0980];
blue   = [0 0.4470 0.7410];
yellow = [0.9290 0.6940 0.1250];

% =========================================================================
% PREPARE INVERSE SIMULATION RELATIONSHIP
% =========================================================================
%
% Original simulation relationship:
%
%       simulation input -> simulated pressure peak
%
% We plot the inverse:
%
%       simulated pressure peak -> required simulation input
%
% This is the relationship evaluated using interp1(...,'pchip').

% Sort simulated peak positions
[sim_peak_sorted, idx] = sort(our_sim);
sim_input_sorted = focus_lengths(idx);

% Remove duplicates if necessary
[sim_peak_unique, unique_idx] = unique(sim_peak_sorted, 'stable');
sim_input_unique = sim_input_sorted(unique_idx);

% =========================================================================
% GENERATE PCHIP CURVE
% =========================================================================

% Interpolation range
z_peak_fine = linspace( ...
    min(sim_peak_unique), ...
    max(sim_peak_unique), ...
    1000);

% Inverse PCHIP interpolation
z_input_pchip = interp1( ...
    sim_peak_unique, ...
    sim_input_unique, ...
    z_peak_fine, ...
    'pchip');

% =========================================================================
% PLOT INVERSE FOCAL-DEPTH CALIBRATION
% =========================================================================

fig = figure;

% Original simulation points
scatter( ...
    our_sim, ...
    focus_lengths, ...
    65, ...
    'filled', ...
    'ColorVariable', blue, ...
    'DisplayName', 'Simulation data');

hold on;

% Hydrophone targets and corresponding corrected simulation inputs
scatter( ...
    tus_device, ...
    corrected_sim_input, ...
    75, ...
    '^', ...
    'LineWidth', 1.5, ...
    'ColorVariable', orange, ...
    'DisplayName', 'Hydrophone targets');

% PCHIP interpolation
plot( ...
    z_peak_fine, ...
    z_input_pchip, ...
    'LineWidth', 1.8, ...
    'Color', yellow, ...
    'DisplayName', 'PCHIP interpolation');



xlabel('Actual focal position (mm)');
ylabel('Simulation focal-depth input (mm)');

title('Inverse focal-depth calibration');

legend('Location', 'northwest');

grid on;
box on;

set(gca, 'FontSize', 12);
set(fig, 'Color', 'w');

% Save figure
if save_outptuts
    saveas(fig, 'inverse_focal_depth_calibration.png');
    savefig(fig, 'inverse_focal_depth_calibration.fig');
end


% =========================================================================
% PRINT CALIBRATION POINTS
% =========================================================================

calibration_table = table( ...
    focus_lengths(:), ...
    tus_device(:), ...
    our_sim(:), ...
    corrected_sim_input(:), ...
    'VariableNames', { ...
    'Commanded_Depth_mm', ...
    'Hydrophone_Peak_mm', ...
    'Uncorrected_Sim_Peak_mm', ...
    'Corrected_Sim_Input_mm'});

fprintf('\n====================================================\n');
fprintf('FOCAL-DEPTH CALIBRATION POINTS\n');
fprintf('====================================================\n\n');

disp(calibration_table);

if save_outputs
    writetable( ...
        calibration_table, ...
        'focal_depth_calibration_points.csv');
end


% =========================================================================
% LOAD VALIDATION RESULTS
% =========================================================================

validation = load('delay_calib_results/delay_calibration_validation.mat');

corrected_input    = validation.corrected_input;
corrected_sim_peak = validation.corrected_sim_peak;

% =========================================================================
% CALCULATE BEFORE / AFTER ERRORS
% =========================================================================

error_before = our_sim - tus_device;
error_after  = corrected_sim_peak - tus_device;

abs_error_before = abs(error_before);
abs_error_after  = abs(error_after);

% =========================================================================
% FINAL THESIS VALIDATION TABLE
% =========================================================================

if save_outputs
    results_table = table( ...
        focus_lengths(:), ...
        tus_device(:), ...
        our_sim(:), ...
        corrected_input(:), ...
        corrected_sim_peak(:), ...
        error_before(:), ...
        error_after(:), ...
        'VariableNames', { ...
        'Commanded_Depth_mm', ...
        'Hydrophone_Peak_mm', ...
        'Uncorrected_Sim_Peak_mm', ...
        'Corrected_Sim_Input_mm', ...
        'Corrected_Sim_Peak_mm', ...
        'Error_Before_mm', ...
        'Error_After_mm'});
    
    fprintf('\n====================================================\n');
    fprintf('DELAY CALIBRATION VALIDATION\n');
    fprintf('====================================================\n\n');
end

disp(results_table);

% =========================================================================
% SUMMARY METRICS
% =========================================================================

MAE_before = mean(abs_error_before);
MAE_after  = mean(abs_error_after);

RMSE_before = sqrt(mean(error_before.^2));
RMSE_after  = sqrt(mean(error_after.^2));

MaxError_before = max(abs_error_before);
MaxError_after  = max(abs_error_after);

fprintf('\nMean absolute error:\n');
fprintf('  Before = %.2f mm\n', MAE_before);
fprintf('  After  = %.2f mm\n', MAE_after);

fprintf('\nRMSE:\n');
fprintf('  Before = %.2f mm\n', RMSE_before);
fprintf('  After  = %.2f mm\n', RMSE_after);

fprintf('\nMaximum absolute error:\n');
fprintf('  Before = %.2f mm\n', MaxError_before);
fprintf('  After  = %.2f mm\n', MaxError_after);

% =========================================================================
% SAVE TABLE
% =========================================================================

writetable( ...
    results_table, ...
    'focal_depth_calibration_validation.csv');

save( ...
    'focal_depth_calibration_results.mat', ...
    'calibration_table', ...
    'results_table', ...
    'MAE_before', ...
    'MAE_after', ...
    'RMSE_before', ...
    'RMSE_after', ...
    'MaxError_before', ...
    'MaxError_after');

fprintf('\nSaved:\n');
fprintf('  inverse_focal_depth_calibration.png\n');
fprintf('  inverse_focal_depth_calibration.fig\n');
fprintf('  focal_depth_calibration_points.csv\n');
fprintf('  focal_depth_calibration_validation.csv\n');
fprintf('  focal_depth_calibration_results.mat\n');