%% Validate delay / focal-depth calibration
% get_validation_delays.m
%
% Re-runs the TRUST simulations using the calibrated focal-depth input 
% (obtained after running delay_calibration.m):
%
%   z_sim_input = polyval(p_delay, z_GUI)
%
% and compares the resulting simulated focal peak against the
% hydrophone-measured focal position.
% 
% Plots Figure 3.4 and Table 3.3 from the Mster's thesis.
%
% author: Maylis Mueller
% date:
% last update:

clear;
clc;
close all;


% =========================================================================
% OUTPUT FOLDERS
% =========================================================================

output_folder_plots = 'delay_calib_results';
if ~exist(output_folder_plots, 'dir')
    mkdir(output_folder_plots);
end

pressure_output_dir = 'delay_calib_results/pressure_fields_corrected';
if ~exist(pressure_output_dir, 'dir')
    mkdir(pressure_output_dir);
end

save_plots = false;
save_data = false;

% =========================================================================
% LOAD DELAY CALIBRATION
% =========================================================================

load('delay_calib_results/delay_calibration.mat', ...
    'focus_lengths', ...
    'tus_device', ...
    'our_sim', ...
    'corrected_sim_input');

num_foci = numel(focus_lengths);

fprintf('\nUsing direct PCHIP focal-depth calibration.\n');


% =========================================================================
% COMPUTATIONAL PARAMETERS
% =========================================================================

t_end = 160e-6;
cfl   = 0.2;


% =========================================================================
% MEDIUM
% =========================================================================

c0   = 1480;
rho0 = 1000;

medium.sound_speed = c0;
medium.density     = rho0;
medium.alpha_coeff = 2.17e-3;
medium.alpha_power = 2;


% =========================================================================
% GRID
% =========================================================================

grid_size_x = 60e-3;
grid_size_y = 60e-3;
grid_size_z = 150e-3;

dx = 0.25e-3;
dy = dx;
dz = dx;

Nx = roundEven(grid_size_x / dx);
Ny = roundEven(grid_size_y / dy);
Nz = roundEven(grid_size_z / dz);

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

kgrid.makeTime(c0, cfl, t_end);

sample_freq = 1 / kgrid.dt;

% Axial coordinate measured from source plane
z_axis_mm = (kgrid.z_vec - kgrid.z_vec(1)) * 1e3;

% Central axis indices
[~, ix0] = min(abs(kgrid.x_vec));
[~, iy0] = min(abs(kgrid.y_vec));


% =========================================================================
% SOURCE PARAMETERS
% =========================================================================

source_f0  = 480e3;
source_amp = 0.1005e6;

row_counts = ...
    [6 10 12 14 16 16 18 18 18 18 18 18 16 16 14 12 10 6];

element_width  = 3.9e-3;
element_height = 3.9e-3;
element_pitch  = 3.08e-3;


% =========================================================================
% SENSOR
% =========================================================================

sensor.mask   = ones(Nx, Ny, Nz);
sensor.record = {'p_max'};


% =========================================================================
% SIMULATION OPTIONS
% =========================================================================

input_args = { ...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off'};


% =========================================================================
% PREALLOCATE
% =========================================================================
corrected_input = zeros(1, num_foci);
corrected_sim_peak = zeros(1, num_foci);
corrected_fwhm_z_mm = nan(1, num_foci);
fwhm_start_mm = nan(1, num_foci);
fwhm_end_mm   = nan(1, num_foci);


% =========================================================================
% RUN CORRECTED SIMULATIONS
% =========================================================================

for i = 1:num_foci

    z_GUI       = focus_lengths(i);
    z_target    = tus_device(i);
    
    z_sim_input = interp1( ...
        focus_lengths, ...
        corrected_sim_input, ...
        z_GUI, ...
        'pchip');

    corrected_input(i) = z_sim_input;

    fprintf('\n============================================\n');
    fprintf('Case %d/%d\n', i, num_foci);
    fprintf('GUI input              = %.2f mm\n', z_GUI);
    fprintf('Hydrophone target      = %.2f mm\n', z_target);
    fprintf('Corrected sim input    = %.2f mm\n', z_sim_input);
    fprintf('============================================\n');

    % Corrected focal point

    focal_point_user = [0 0 z_sim_input] * 1e-3;

    focal_point_grid = [ ...
        focal_point_user(1), ...
        focal_point_user(2), ...
        kgrid.z_vec(1) + focal_point_user(3)];

    % Create source

    [source, signal_info] = create_source( ...
        focal_point_grid, ...
        row_counts, ...
        element_height, ...
        element_width, ...
        element_pitch, ...
        source_amp, ...
        sample_freq, ...
        source_f0, ...
        c0, ...
        kgrid, ...
        'NumCycles', 30);

    % Run simulation

    sensor_data = kspaceFirstOrder3DG( ...
        kgrid, ...
        medium, ...
        source, ...
        sensor, ...
        input_args{:}, ...
        'DataCast', 'single', ...
        'PlotScale', [-1 1] * source_amp, ...
        'PlotSim', false);

    % Reshape pressure field

    p_max = reshape(sensor_data.p_max, [Nx Ny Nz]);

    % Save full pressure field for this focal depth

    pressure_filename = fullfile( ...
        pressure_output_dir, ...
        sprintf('pressure_field_GUI_%03dmm.mat', round(z_GUI)));
    
    if save_data
        save( ...
            pressure_filename, ...
            'p_max', ...
            'z_GUI', ...
            'z_target', ...
            'z_sim_input', ...
            'dx', 'dy', 'dz', ...
            'source_amp', ...
            'source_f0', ...
            '-v7.3');
        
        fprintf('Pressure field saved: %s\n', pressure_filename);
    end

    % Extract on-axis axial profile
   
    axial_profile = double(squeeze(p_max(ix0, iy0, :)));
    
    % Search around hydrophone target
    
    search_width = 30;     % [mm]
    
    z_min = max(0, z_target - search_width);
    z_max = min(grid_size_z * 1e3, z_target + search_width);
    
    search_mask = ...
        z_axis_mm >= z_min & ...
        z_axis_mm <= z_max;
    
    search_ind = find(search_mask);
    
    z_search = z_axis_mm(search_ind);
    p_search = axial_profile(search_ind);
    
    % Find corrected simulated focal peak
    
    [~, idx_peak_local] = max(p_search);
    
    % Convert index in search region back to full axial-profile index
    peak_ind = search_ind(idx_peak_local);
    
    corrected_sim_peak(i) = z_axis_mm(peak_ind);
    
    % Axial FWHM around corrected focal peak
    
    [corrected_fwhm_z_mm(i), ...
     fwhm_start_mm(i), ...
     fwhm_end_mm(i)] = ...
        get_fwhm( ...
            z_axis_mm, ...
            axial_profile, ...
            peak_ind);
    
    fprintf('Corrected sim peak     = %.2f mm\n', ...
        corrected_sim_peak(i));
    
    fprintf('Error vs hydrophone    = %+.2f mm\n', ...
        corrected_sim_peak(i) - z_target);
    
    fprintf('FWHM_z                 = %.2f mm\n', ...
        corrected_fwhm_z_mm(i));
    
    fprintf('FWHM limits            = %.2f to %.2f mm\n', ...
        fwhm_start_mm(i), ...
        fwhm_end_mm(i));
end


% =========================================================================
% COMPUTE ERRORS
% =========================================================================

uncorrected_error = our_sim - tus_device;
corrected_error   = corrected_sim_peak - tus_device;

abs_uncorrected_error = abs(uncorrected_error);
abs_corrected_error   = abs(corrected_error);


% =========================================================================
% RESULTS TABLE
% =========================================================================

results_validation = table( ...
    focus_lengths(:), ...
    tus_device(:), ...
    our_sim(:), ...
    uncorrected_error(:), ...
    corrected_input(:), ...
    corrected_sim_peak(:), ...
    corrected_error(:), ...
    abs_corrected_error(:), ...
    fwhm_start_mm(:), ...
    fwhm_end_mm(:), ...
    corrected_fwhm_z_mm(:), ...
    'VariableNames', { ...
    'GUI_Input_mm', ...
    'Hydrophone_Focus_mm', ...
    'Uncorrected_Sim_Peak_mm', ...
    'Uncorrected_Error_mm', ...
    'Corrected_Sim_Input_mm', ...
    'Corrected_Sim_Peak_mm', ...
    'Corrected_Error_mm', ...
    'Abs_Corrected_Error_mm', ...
    'FWHM_Start_mm', ...
    'FWHM_End_mm', ...
    'Corrected_FWHM_Z_mm'});

disp(results_validation);


% =========================================================================
% SUMMARY METRICS
% =========================================================================

mean_abs_error_before = mean(abs_uncorrected_error);
mean_abs_error_after  = mean(abs_corrected_error);

max_abs_error_before = max(abs_uncorrected_error);
max_abs_error_after  = max(abs_corrected_error);

rmse_before = sqrt(mean(uncorrected_error.^2));
rmse_after  = sqrt(mean(corrected_error.^2));

fprintf('\n====================================================\n');
fprintf('DELAY CALIBRATION VALIDATION\n');
fprintf('====================================================\n');

fprintf('Mean absolute error before: %.3f mm\n', ...
    mean_abs_error_before);

fprintf('Mean absolute error after : %.3f mm\n', ...
    mean_abs_error_after);

fprintf('\n');

fprintf('Maximum absolute error before: %.3f mm\n', ...
    max_abs_error_before);

fprintf('Maximum absolute error after : %.3f mm\n', ...
    max_abs_error_after);

fprintf('\n');

fprintf('RMSE before: %.3f mm\n', rmse_before);
fprintf('RMSE after : %.3f mm\n', rmse_after);

fprintf('====================================================\n');


% =========================================================================
% PLOT
% =========================================================================

% colors
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
light_blue = [0.5843 0.8157 0.9882];

x = results_validation.GUI_Input_mm;
hyd = results_validation.Hydrophone_Focus_mm;
sim = results_validation.Uncorrected_Sim_Peak_mm;
sim_corr = results_validation.Corrected_Sim_Peak_mm;

fig = figure; 
hold on;

lims = [min(x) max(x)];

plot(lims, lims, 'k:', 'LineWidth', 1.5, ...
    'DisplayName', 'Commanded = realised');

plot(x, hyd, '^-', 'Color', orange, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'Hydrophone');

plot(x, sim, 'o-', 'Color', blue, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'Uncorrected simulation');

plot(x, sim_corr, 's-', 'Color', light_blue, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'Corrected simulation');

xlabel('Commanded focal depth [mm]', 'FontSize',14);
ylabel('Realised focal depth [mm]', 'FontSize',14);
legend('Location', 'southeast', 'FontSize',14);
grid on; box on;

% Same scale on x and y makes deviation from x = y easier to interpret
xlim(lims);
ylim(lims);
axis square;

if save_plots
    saveas(fig, fullfile(output_folder_plots, 'delay_calibration_validation.png'));
    savefig(fig, fullfile(output_folder_plots, 'delay_calibration_validation.fig'));
end


% =========================================================================
% SAVE
% =========================================================================

if save_data
    save( ...
        fullfile(output_folder_plots, 'delay_calibration_validation.mat'), ...
        'results_validation', ...
        'focus_lengths', ...
        'tus_device', ...
        'corrected_input', ...
        'corrected_sim_peak', ...
        'uncorrected_error', ...
        'corrected_error', ...
        'corrected_fwhm_z_mm', ...
        'fwhm_start_mm', ...
        'fwhm_end_mm', ...
        'mean_abs_error_before', ...
        'mean_abs_error_after', ...
        'max_abs_error_before', ...
        'max_abs_error_after', ...
        'rmse_before', ...
        'rmse_after');
    
    writetable( ...
        results_validation, ...
        'delay_calibration_validation.csv');
    
    fprintf('\nResults saved to:\n');
    fprintf('  delay_calibration_validation.mat\n');
    fprintf('  delay_calibration_validation.csv\n');
end