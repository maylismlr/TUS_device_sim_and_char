%% Delay / focal-depth calibration
% delay_calibration.m
% Determines the corrected focal-depth input required in the simulation
% so that the simulated axial pressure maximum matches the focal position
% measured with the TRUST device.
%
% Calibration logic:
%
%   TRUST GUI input
%       -> measured TRUST focal position
%       -> inverse of simulated input/actual-focus relationship
%       -> corrected simulation input
%
% author: Maylis Mueller
% date:
% last update:

clear;
clc;
close all;

output_dir = 'delay_calib_results';

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% =========================================================================
% USER-DEFINED PARAMETERS
% =========================================================================

% Computational parameters
t_end = 160e-6;        % total simulation time [s]
cfl   = 0.2;           % CFL number

% Medium: water
c0   = 1480;           % sound speed [m/s]
rho0 = 1000;           % density [kg/m^3]

medium.sound_speed = c0;
medium.density     = rho0;
medium.alpha_coeff = 2.17e-3;
medium.alpha_power = 2;

% Grid
grid_size_x = 60e-3;   % [m]
grid_size_y = 60e-3;   % [m]
grid_size_z = 150e-3;  % [m]

dx = 0.25e-3;
dy = dx;
dz = dx;

save_plots = false;
save_data = false;


% =========================================================================
% CREATE GRID
% =========================================================================

Nx = roundEven(grid_size_x / dx);
Ny = roundEven(grid_size_y / dy);
Nz = roundEven(grid_size_z / dz);

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% Time array
kgrid.makeTime(c0, cfl, t_end);

sample_freq = 1 / kgrid.dt;

% Axial coordinate expressed as distance from the source plane
z_axis_mm = (kgrid.z_vec - kgrid.z_vec(1)) * 1e3;

% Central x/y grid indices
[~, ix0] = min(abs(kgrid.x_vec));
[~, iy0] = min(abs(kgrid.y_vec));


% =========================================================================
% SOURCE PARAMETERS
% =========================================================================

source_f0     = 480e3;       % [Hz]
source_amp    = 0.1005e6;    % [Pa], arbitrary for focal-position calibration
source_cycles = 30;

% TRUST array geometry
row_counts = ...
    [6 10 12 14 16 16 18 18 18 18 18 18 16 16 14 12 10 6];

% Calibrated effective element size
element_width  = 3.9e-3;     % [m]
element_height = 3.9e-3;     % [m]

% Physical element pitch
element_pitch = 3.08e-3;     % [m]

% Temporal parameters
PRF_period = 15e-6;
num_bursts = 1;


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
% REFERENCE DATA
% =========================================================================

% Focal depths requested on the TRUST device
% These must correspond one-to-one with the hydrophone profiles.
focus_lengths = ...
    [20 30 40 50 60 70 80 90 100 110 120];

% Focal positions measured from the corresponding hydrophone profiles
tus_device = ...
    [19 30 40 51 59 69 75 83 89 97 105];

assert(numel(focus_lengths) == numel(tus_device), ...
    'focus_lengths and tus_device must contain the same number of points.');

num_foci = numel(focus_lengths);


% =========================================================================
% RUN UNCORRECTED SIMULATIONS
% =========================================================================

% Actual focal positions produced by the simulation
our_sim = zeros(1, num_foci);

for i = 1:num_foci

    z_input = focus_lengths(i);

    fprintf('\n========================================\n');
    fprintf('Focus %d/%d\n', i, num_foci);
    fprintf('Requested simulation focus = %.1f mm\n', z_input);
    fprintf('========================================\n');

    % Requested focal point

    focal_point_user = [0 0 z_input] * 1e-3;

    % Convert depth measured from source plane to k-Wave coordinates
    focal_point_grid = [ ...
        focal_point_user(1), ...
        focal_point_user(2), ...
        kgrid.z_vec(1) + focal_point_user(3)];

    % Create source

    [source, signal_info] = create_source(focal_point_grid, row_counts, ...
            element_height, element_width, element_pitch, ...
            source_amp, sample_freq, source_f0, c0, kgrid, ...
            'NumCycles', 30);

    % Run simulation

    sensor_data = kspaceFirstOrder3DG( ...
        kgrid, medium, source, sensor, ...
        input_args{:}, ...
        'DataCast', 'single', ...
        'PlotScale', [-1 1] * source_amp, ...
        'PlotSim', false);

    % Reshape pressure field

    p_max = reshape(sensor_data.p_max, [Nx Ny Nz]);

    % Extract on-axis pressure profile

    axial_profile = squeeze(p_max(ix0, iy0, :));

    % Search for axial maximum around expected focus

    % Prevent the source region or unrelated maxima from being selected.
    search_width = 30;       % [mm]

    z_min = max(0, z_input - search_width);
    z_max = min(grid_size_z * 1e3, z_input + search_width);

    search_mask = ...
        z_axis_mm >= z_min & ...
        z_axis_mm <= z_max;

    z_search = z_axis_mm(search_mask);
    p_search = axial_profile(search_mask);

    % Find simulated focal position directly on the k-Wave grid
    [~, idx_peak] = max(p_search);
    
    our_sim(i) = z_search(idx_peak);

    fprintf('Actual simulated focus = %.2f mm\n', our_sim(i));
    fprintf('Axial shift = %.2f mm\n', our_sim(i) - z_input);

end


% =========================================================================
% DISPLAY UNCORRECTED RESULTS
% =========================================================================

results_uncorrected = table( ...
    focus_lengths(:), ...
    tus_device(:), ...
    our_sim(:), ...
    our_sim(:) - focus_lengths(:), ...
    our_sim(:) - tus_device(:), ...
    'VariableNames', { ...
    'GUI_Input_mm', ...
    'Device_Focus_mm', ...
    'Sim_Focus_mm', ...
    'Sim_Shift_From_Input_mm', ...
    'Sim_Error_vs_Device_mm'});

disp(results_uncorrected);


% =========================================================================
% PLOT 1
% SIMULATION INPUT -> ACTUAL SIMULATED FOCUS
% =========================================================================

fig1 = figure;

scatter(focus_lengths, our_sim, 50, 'filled');
hold on;

% Identity line
limits = [ ...
    min([focus_lengths our_sim]), ...
    max([focus_lengths our_sim])];

plot(limits, limits, '--', 'LineWidth', 1.2);

xlabel('Simulation focal-depth input [mm]');
ylabel('Actual simulated focal position [mm]');
title('Uncorrected simulation focal-depth behaviour');

legend( ...
    'Simulation', ...
    'Input = actual focus', ...
    'Location', 'best');

grid on;
box on;

if save_plots
    saveas(fig1, fullfile(output_dir, 'simulation_input_vs_actual_focus.png'));
    savefig(fig1, fullfile(output_dir, 'simulation_input_vs_actual_focus.fig'));
end


% =========================================================================
% INVERT SIMULATION BEHAVIOUR
% =========================================================================
%
% We know:
%
%   focus_lengths -> our_sim
%
% i.e.:
%
%   simulation input -> actual simulated focus
%
% We now invert this relationship:
%
%   desired actual focus -> required simulation input
%
% The desired actual focus is the hydrophone-measured TRUST focus.

% Sort simulated actual focal positions before using interp1
[our_sim_sorted, sort_idx] = sort(our_sim);

focus_lengths_sorted = focus_lengths(sort_idx);

% Remove duplicate simulated focal positions if any
[our_sim_unique, unique_idx] = unique( ...
    our_sim_sorted, ...
    'stable');

focus_lengths_unique = focus_lengths_sorted(unique_idx);

% Check whether hydrophone values require extrapolation
if min(tus_device) < min(our_sim_unique) || ...
   max(tus_device) > max(our_sim_unique)

    warning([ ...
        'Some TRUST focal positions lie outside the range covered by ' ...
        'the simulated actual focal positions. Extrapolation will be used.']);

end

% Required simulation input to reproduce each hydrophone focal position
corrected_sim_input = interp1( ...
    our_sim_unique, ...
    focus_lengths_unique, ...
    tus_device, ...
    'pchip', ...
    'extrap');


% =========================================================================
% FINAL GUI -> CORRECTED SIMULATION INPUT CALIBRATION
% =========================================================================
%
% At this stage:
%
% focus_lengths        = TRUST GUI/requested focus
% corrected_sim_input  = input our simulation should actually receive
%
% Fit a linear equation for convenient implementation in the model:
%
% z_sim,input = a * z_GUI + b

p_delay = polyfit( ...
    focus_lengths, ...
    corrected_sim_input, ...
    1);

delay_slope     = p_delay(1);
delay_intercept = p_delay(2);

fprintf('\n\n');
fprintf('====================================================\n');
fprintf('FINAL DELAY / FOCAL-DEPTH CALIBRATION\n');
fprintf('====================================================\n');

fprintf( ...
    'z_sim_input = %.6f * z_GUI %+.6f mm\n', ...
    delay_slope, ...
    delay_intercept);

fprintf('====================================================\n');


% =========================================================================
% PLOT 2
% CORRECTED SIMULATION INPUT
% =========================================================================

x_fit = linspace( ...
    min(focus_lengths), ...
    max(focus_lengths), ...
    500);

y_fit = polyval(p_delay, x_fit);

fig2 = figure;

scatter( ...
    focus_lengths, ...
    corrected_sim_input, ...
    50, ...
    'filled');

hold on;

plot( ...
    x_fit, ...
    y_fit, ...
    'LineWidth', 1.5);

% No-correction reference
plot( ...
    x_fit, ...
    x_fit, ...
    '--', ...
    'LineWidth', 1.2);

xlabel('TRUST requested focal depth [mm]');
ylabel('Corrected simulation focal-depth input [mm]');
title('Focal-depth calibration');

legend( ...
    'Calibration points', ...
    'Linear calibration', ...
    'No correction', ...
    'Location', 'best');

grid on;
box on;

eqn_str = sprintf( ...
    'z_{sim,input} = %.4f z_{GUI} %+.4f', ...
    delay_slope, ...
    delay_intercept);

text( ...
    0.05, ...
    0.92, ...
    eqn_str, ...
    'Units', 'normalized', ...
    'FontSize', 11);

if save_plots
    saveas(fig2, fullfile(output_dir, 'final_delay_calibration.png'));
    savefig(fig2, fullfile(output_dir, 'final_delay_calibration.fig'));
end


% =========================================================================
% CALIBRATION TABLE
% =========================================================================

correction_mm = corrected_sim_input - focus_lengths;

results_calibration = table( ...
    focus_lengths(:), ...
    tus_device(:), ...
    our_sim(:), ...
    corrected_sim_input(:), ...
    correction_mm(:), ...
    'VariableNames', { ...
    'GUI_Input_mm', ...
    'Device_Focus_mm', ...
    'Uncorrected_Sim_Focus_mm', ...
    'Corrected_Sim_Input_mm', ...
    'Input_Correction_mm'});

disp(results_calibration);


% =========================================================================
% PLOT 3
% DEVICE FOCUS VS UNCORRECTED SIMULATION FOCUS
% =========================================================================

fig3 = figure;

scatter( ...
    tus_device, ...
    our_sim, ...
    50, ...
    'filled');

hold on;

limits = [ ...
    min([tus_device our_sim]), ...
    max([tus_device our_sim])];

plot( ...
    limits, ...
    limits, ...
    '--', ...
    'LineWidth', 1.2);

xlabel('Measured TRUST focal position [mm]');
ylabel('Uncorrected simulated focal position [mm]');
title('Measured vs uncorrected simulated focal position');

legend( ...
    'Calibration cases', ...
    'Perfect agreement', ...
    'Location', 'best');

grid on;
box on;

if save_plots
    saveas(fig3, fullfile(output_dir, 'device_vs_uncorrected_sim_focus.png'));
    savefig(fig3, fullfile(output_dir, 'device_vs_uncorrected_sim_focus.fig'));
end


% =========================================================================
% OPTIONAL LINEAR CHARACTERISATION OF UNCORRECTED ERROR
% =========================================================================
%
% This is useful for reporting/visualisation, but is NOT what is used
% to obtain the corrected simulation inputs above.

p_device_sim = polyfit( ...
    our_sim, ...
    tus_device, ...
    1);

fprintf('\nUncorrected simulated-focus -> device-focus relation:\n');

fprintf( ...
    'z_device = %.6f * z_sim_actual %+.6f mm\n', ...
    p_device_sim(1), ...
    p_device_sim(2));


% =========================================================================
% SAVE CALIBRATION
% =========================================================================

if save_data
    save( ...
        fullfile(output_dir, 'delay_calibration.mat'), ...
        'focus_lengths', ...
        'tus_device', ...
        'our_sim', ...
        'corrected_sim_input', ...
        'correction_mm', ...
        'p_delay', ...
        'delay_slope', ...
        'delay_intercept', ...
        'p_device_sim', ...
        'results_uncorrected', ...
        'results_calibration');
    
    fprintf('\nCalibration saved to delay_calibration.mat\n');
end