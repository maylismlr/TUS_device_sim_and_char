%% Validate input pressure calibration
% plot_input_p0_validation.m
%
% Re-runs the TRUST simulations using the calibrated focal-depth input (initialized by
% calling the trust_reference_data() function and compares the resulting 
% simulated focal peak against the hydrophone-measured focal position.
%
% Plots Figure 3.6 from Master's thesis.
%
% author: Maylis Mueller
% date:
% last update:

clear;
clc;
close all;

addpath '/data/tu_muellerm/utils/'
addpath '/data/hu_muellerm/MATLAB Add-Ons/Collections/k-Wave/k-Wave'

% for input pressure calibration
p_ref = trust_reference_data();

save_data=false;

% =========================================================================
% PRESSURE FIELD OUTPUT FOLDER
% =========================================================================

pressure_output_dir = 'input_calib_pressure_fields';

if ~exist(pressure_output_dir, 'dir')
    mkdir(pressure_output_dir);
end

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
source_volt = 30;

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
focus_lengths = p_ref.z_hyd;
num_foci = numel(focus_lengths);

corrected_input_p0 = zeros(1, num_foci);
target_peff = zeros(1, num_foci);
corrected_pmax = zeros(1, num_foci);
pressure_error = zeros(1, num_foci);
pressure_epsilon = zeros(1, num_foci);


% =========================================================================
% RUN CORRECTED SIMULATIONS
% =========================================================================

for i = 1:num_foci

    z_GUI       = focus_lengths(i);
    
    fprintf('\n============================================\n');
    fprintf('Case %d/%d\n', i, num_foci);
    fprintf('GUI input              = %.2f mm\n', z_GUI);
    fprintf('============================================\n');

    % Corrected focal point
    focal_point_user = [0 0 z_GUI] * 1e-3;

    focal_point_grid = [ ...
        focal_point_user(1), ...
        focal_point_user(2), ...
        kgrid.z_vec(1) + focal_point_user(3)];

    % Create source
    % adding input pressure calib
    [source_amp, p0_info] = trust_p0(source_volt, z_GUI, p_ref);

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

    % Isolate quantities of interest
    p0_target = p0_info.p_target;

    % Reshape pressure field
    p_max = reshape(sensor_data.p_max, [Nx Ny Nz]);

    if save_data
        % Save full pressure field for this focal depth
        pressure_filename = fullfile( ...
            pressure_output_dir, ...
            sprintf('pressure_field_GUI_%03dmm.mat', round(z_GUI)));
        
        save( ...
            pressure_filename, ...
            'p_max', ...
            'z_GUI', ...
            'p0_target', ...
            'source_amp', ...
            'dx', 'dy', 'dz', ...
            'source_f0', ...
            '-v7.3');
        
        fprintf('Pressure field saved: %s\n', pressure_filename);
    end

    % Extract on-axis axial profile
    axial_profile = double(squeeze(p_max(ix0, iy0, :)));
    
    % Search around hydrophone target
    search_width = 30;     % [mm]
    
    z_min = max(0, z_GUI - search_width);
    z_max = min(grid_size_z * 1e3, z_GUI + search_width);
    
    search_mask = ...
        z_axis_mm >= z_min & ...
        z_axis_mm <= z_max;
    
    search_ind = find(search_mask);
    
    z_search = z_axis_mm(search_ind);
    p_search = axial_profile(search_ind);
    
    % Find corrected simulated max pressure
    [pmax, idx_peak_local] = max(p_search);

    corrected_input_p0(i) = source_amp;
    corrected_pmax(i) = pmax;
    
    target_peff(i) = p0_info.p_target;
    
    pressure_error(i) = corrected_pmax(i) - target_peff(i);
    pressure_epsilon(i) = 100 * pressure_error(i) / target_peff(i);
    
    fprintf('Input source pressure = %.3f kPa\n', source_amp / 1e3);
    fprintf('Target focal pressure = %.3f MPa\n', target_peff(i) / 1e6);
    fprintf('Simulated focal p_max = %.3f MPa\n', corrected_pmax(i) / 1e6);
    fprintf('Pressure error        = %+.3f MPa (%+.2f%%)\n', ...
        pressure_error(i) / 1e6, pressure_epsilon(i));
   
end



% =========================================================================
% RESULTS TABLE
% =========================================================================

results_validation = table( ...
    focus_lengths(:), ...
    target_peff(:), ...
    corrected_input_p0(:), ...
    corrected_pmax(:), ...
    pressure_error(:), ...
    pressure_epsilon(:), ...
    'VariableNames', { ...
    'GUI_Input_mm', ...
    'Hydrophone_P_eff', ...
    'Corrected_p0_Input', ...
    'Corrected_Sim_pmax', ...
    'Corrected_Error', ...
    'Error_epsilon'});

disp(results_validation);


% =========================================================================
% PLOT
% =========================================================================
blue   = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];

fig = figure;
hold on;

plot(focus_lengths, target_peff / 1e6, '^-', ...
    'Color', orange, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6, ...
    'DisplayName', 'Hydrophone target');

plot(focus_lengths, corrected_pmax / 1e6, 'o-', ...
    'Color', blue, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6, ...
    'DisplayName', 'Calibrated simulation');

xlabel('Commanded focal depth [mm]', 'FontSize',14);
ylabel('Focal pressure [MPa]', 'FontSize',14);
legend('Location', 'best', 'FontSize',14);
grid on;
box on;

savefig(fig, 'pressure_calibration_validation.fig');
exportgraphics(fig, 'pressure_calibration_validation.png', 'Resolution', 300);

% =========================================================================
% SAVE
% =========================================================================

pressure_results = table( ...
    focus_lengths(:), ...
    corrected_input_p0(:) / 1e3, ...
    target_peff(:) / 1e6, ...
    corrected_pmax(:) / 1e6, ...
    pressure_error(:) / 1e6, ...
    pressure_epsilon(:), ...
    'VariableNames', { ...
    'CommandedDepth_mm', ...
    'InputPressure_kPa', ...
    'TargetPressure_MPa', ...
    'SimulatedPressure_MPa', ...
    'PressureError_MPa', ...
    'PressureEpsilon'});

disp(pressure_results);

writetable(pressure_results, 'pressure_calibration_validation.csv');

fprintf('\nPRESSURE CALIBRATION VALIDATION\n');
fprintf('Mean absolute error = %.3f MPa\n', mean(abs(pressure_error)) / 1e6);
fprintf('Mean absolute percentage error = %.2f %%\n', mean(abs(pressure_epsilon)));
fprintf('Maximum absolute percentage error = %.2f %%\n', max(abs(pressure_epsilon)));
fprintf('RMSE = %.3f MPa\n', sqrt(mean(pressure_error.^2)) / 1e6);