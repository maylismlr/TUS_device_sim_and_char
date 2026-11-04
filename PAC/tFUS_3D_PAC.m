%% TRUST 3-D pCT simulation with phase-aberration correction
% Based on tFUS_3D_MR.m.
% Runs the whole PAC pipeline. To check if it worked, replaced medium_water
% by medium_skull in the estimate_pac_phase_reversal call. Gave out same
% simulations as without PAC and flat progress plots.
%
% Author: Maylis Mueller

clearvars;
addpath '/data/tu_muellerm/utils/'
addpath '/data/hu_muellerm/MATLAB Add-Ons/Collections/k-Wave/k-Wave'

% =========================================================================
% USER SETTINGS
% =========================================================================
% muCT
name = "muCT";
a = "";

% MR_HR_a2862
% name = "MR_HR";
% a = "_a2862";

path_medium = sprintf('/data/tu_muellerm_cloud/skulls_FUN/skulls_params/final/%s.nrrd%s_0.25mm_params.mat', name, a);

path_delay_calibration = ...
    '/data/tu_muellerm/utils/delay_calibration.mat'; % adjust if needed
add_delays = false;


target_point_user = [0, 0, 45] * 1e-3;

x_steer_mm = target_point_user(1) * 1e3;
y_steer_mm = target_point_user(2) * 1e3;
z_focus_mm = target_point_user(3) * 1e3;

% PAC estimation consists of one target-source run in water and one with the
% skull. Set false after pac_file has been generated for this target/skull.
run_pac_estimation = true;
apply_pac          = true;  % false gives the otherwise identical no-PAC case

pac_file = sprintf( ...
    '/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/pac/pac_%s%s_%gmm_%gmm_focus_reciprocal.mat', name, a, x_steer_mm, z_focus_mm);

output_uncorrected = sprintf( ...
    '/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/%s%s_p_max_%gmm_%gmm_noPAC.mat', name, a, x_steer_mm, z_focus_mm);
output_corrected = sprintf( ...
    '/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/%s%s_p_max_%gmm_%gmm_PAC.mat', name, a, x_steer_mm, z_focus_mm);

save_data=false;
save_plots=false;

t_end = 120e-6;
cfl = 0.2;
c0 = 1480;
rho0 = 1000;

dx = 0.25e-3;
dy = dx;
dz = dx;
grid_size_x = 56e-3;
grid_size_y = 56e-3;
grid_size_z = 70e-3;

Nx = roundEven(grid_size_x / dx);
Ny = roundEven(grid_size_y / dy);
Nz = roundEven(grid_size_z / dz);


source_f0 = 4.8e5;
source_volt = 30;
source_cycles = 30;

row_counts = [6 10 12 14 16 16 18 18 18 18 18 18 16 16 14 12 10 6];

element_width = 3.9e-3;
element_height = element_width;
element_pitch = 3.08e-3;

% =========================================================================
% GRID AND SKULL MEDIUM
% =========================================================================

% saved elsewhere because too heavy !
tmp = load(path_medium);
seg = tmp.seg();

alpha_water = 2.17e-3;      % dB/(MHz^2 cm)
mask_water = seg(:,:,1:Nz) == 0;


if contains(name, "MR")
    medium_skull.sound_speed = tmp.speed_of_sound(:,:,1:Nz);
    medium_skull.density = tmp.density(:,:,1:Nz);
    medium_skull.alpha_coeff = tmp.alpha_coeff(:,:,1:Nz);
    medium_skull.alpha_coeff(mask_water) = alpha_water;
    medium_skull.alpha_power = 2;
elseif contains(name, "muCT")
    medium_skull.sound_speed = ones(Nx,Ny,Nz);
    seg = seg(1:Nx, 1:Ny, 1:Nz);
    medium_skull.sound_speed(seg==0) = 1480.0;
    medium_skull.sound_speed(seg==1) = 3100.0;
    medium_skull.density = ones(Nx,Ny,Nz);
    medium_skull.density(seg==0) = 1000.0;
    medium_skull.density(seg==1) = 2000.0;
    medium_skull.alpha_coeff = ones(Nx,Ny,Nz);
    medium_skull.alpha_coeff(seg==0) = 2.17e-3;      % dB/(MHz^2 cm)
    medium_skull.alpha_coeff(seg==1) = 2.7;
    medium_skull.alpha_power = 2;
else
    error("Enter a valid skull name")
end


% The reference must be on the identical grid and use the same water model.
medium_water.sound_speed = c0 * ones(Nx, Ny, Nz);
medium_water.density = rho0 * ones(Nx, Ny, Nz);
medium_water.alpha_coeff = alpha_water * ones(Nx, Ny, Nz);
medium_water.alpha_power = 2;

c_ref = max(medium_skull.sound_speed(:));   % same for both

medium_skull.sound_speed_ref = c_ref;
medium_water.sound_speed_ref = c_ref;

if size(seg, 1) ~= Nx || size(seg, 2) ~= Ny || size(seg, 3) < Nz
    error(['The segmentation must already have size Nx-by-Ny and contain ' ...
        'at least Nz axial slices.']);
end

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
kgrid.makeTime(c_ref, cfl, t_end);
sample_freq = 1 / kgrid.dt;

% Desired physical target: this becomes the virtual source location during
% the two reciprocal PAC simulations.
target_point_grid = [target_point_user(1), target_point_user(2), ...
    kgrid.z_vec(1) + target_point_user(3)];

% =========================================================================
% FREE-FIELD DEPTH CALIBRATION
% =========================================================================

z_GUI_mm = target_point_user(3) * 1e3;

% Start from the requested target.
source_focus_user = target_point_user;

if add_delays

    delay_cal = load(path_delay_calibration, ...
        'focus_lengths', 'corrected_sim_input');

    % Check that requested depth is within calibrated range
    if z_GUI_mm < min(delay_cal.focus_lengths) || ...
            z_GUI_mm > max(delay_cal.focus_lengths)

        error(['Requested focal depth %.1f mm is outside the calibrated ' ...
            'range [%.1f, %.1f] mm.'], ...
            z_GUI_mm, ...
            min(delay_cal.focus_lengths), ...
            max(delay_cal.focus_lengths));
    end

    % Convert requested GUI depth to corrected simulation input
    z_sim_input_mm = interp1( ...
        delay_cal.focus_lengths, ...
        delay_cal.corrected_sim_input, ...
        z_GUI_mm, ...
        'pchip');

    % Only z is corrected
    source_focus_user(3) = z_sim_input_mm * 1e-3;

end

% Convert source focus to k-Wave coordinates
source_focus_grid = [source_focus_user(1), source_focus_user(2), ...
    kgrid.z_vec(1) + source_focus_user(3)];

% Amplitude calibration remains indexed by the requested GUI depth.
p_ref = trust_reference_data();
[source_amp, p0_info] = trust_p0(source_volt, z_GUI_mm, p_ref);


solver_args = { ...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off', ...
    'DataCast', 'single'};

% =========================================================================
% ESTIMATE OR LOAD THE ELEMENT-WISE PAC (TWO RECIPROCAL SIMULATIONS)
% =========================================================================

if run_pac_estimation
    pac_amp = 1e5; % virtual target-source amplitude [Pa], 
            % technically could be source_amp but so the distinction is clear

    pac = estimate_pac_phase_reversal( ...
        kgrid, medium_water, medium_skull, target_point_grid, ...
        row_counts, element_height, element_width, element_pitch, ...
        source_f0, pac_amp, solver_args, ...
        'NumCycles', source_cycles, ...
        'SaveFile', pac_file, ...
        'PlotProgress', true, 'SaveProgressPlot', false);
else
    tmp_pac = load(pac_file, 'pac');
    pac = tmp_pac.pac;
end

if apply_pac
    element_pac_delay = pac.correction_delay;
    output_file = output_corrected;
else
    element_pac_delay = zeros(1, sum(row_counts));
    output_file = output_uncorrected;
end

% =========================================================================
% CREATE SOURCE AND RUN THE FULL SKULL SIMULATION
% =========================================================================

[source, signal_info, applied_delays] = create_source( ...
    source_focus_grid, row_counts, element_height, element_width, ...
    element_pitch, source_amp, sample_freq, source_f0, c0, kgrid, ...
    'PACDelay', element_pac_delay, 'NumCycles', source_cycles);


sensor.mask = true(Nx, Ny, Nz);
sensor.record = {'p_max'};

sensor_data = kspaceFirstOrder3DG( ...
    kgrid, medium_skull, source, sensor, solver_args{:}, ...
    'PlotScale', [-1, 1] * source_amp);

p_max = reshape(sensor_data.p_max, [Nx, Ny, Nz]);

if save_data
    save(output_file, 'p_max', 'pac', 'applied_delays', ...
        'target_point_user', 'target_point_grid', ...
        'source_focus_user', 'source_focus_grid', ...
        'source_amp', 'p0_info', 'z_GUI_mm', 'add_delays', '-v7.3');
end

% =========================================================================
% VISUALISATION
% =========================================================================

[~, y_target_index] = min(abs(kgrid.y_vec - target_point_grid(2)));
prop_vec_shifted = 1e3 * (kgrid.z_vec - kgrid.z_vec(1));

fig_field = figure;
imagesc(1e3 * kgrid.x_vec, prop_vec_shifted, ...
    (1e-6 * squeeze(p_max(:, y_target_index, :))).');
xlabel('x-position [mm]');
ylabel('z-position [mm]');
axis image;
title(sprintf('%s pressure field, PAC = %d', name, apply_pac));
cb = colorbar;
title(cb, '[MPa]');

if save_plots

    fig_name = sprintf('pressure_field_pac_%d_steering_%gmm_%gmm_focus', apply_pac, x_steer_mm, z_focus_mm);

    savefig(fig_field, fullfile([fig_name '.fig']));
    exportgraphics(fig_field, fullfile([fig_name '.png']), 'Resolution', 300);

end


% =========================================================================
% PLOT ELEMENT DELAYS
% =========================================================================

positions_mm = pac.positions(:,1:2) * 1e3;

geometric_delay_us = applied_delays.geometric(:) * 1e6;
corrected_delay_us = applied_delays.total(:) * 1e6;

element_x_mm = element_height * 1e3;
element_y_mm = element_width * 1e3;

delay_max = max([geometric_delay_us; corrected_delay_us]);

fig_delays = figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% A: GEOMETRIC DELAYS

ax1 = nexttile;
hold on;

for n = 1:size(positions_mm,1)
    x = positions_mm(n,1) + element_x_mm/2 * [-1 1 1 -1];
    y = positions_mm(n,2) + element_y_mm/2 * [-1 -1 1 1];
    patch(x, y, geometric_delay_us(n), 'EdgeColor', 'w');
end

axis equal;
xlabel('x [mm]', 'FontSize', 14);
ylabel('y [mm]', 'FontSize', 14);
title('Geometrical delays', 'FontSize', 17);
clim([0 delay_max]);
cb = colorbar;
cb.FontSize = 14;
ylabel(cb, 'Delay [\mus]', 'FontSize', 14);
grid on;

% B: ABERRATION-CORRECTED DELAYS

ax2 = nexttile;
hold on;

for n = 1:size(positions_mm,1)
    x = positions_mm(n,1) + element_x_mm/2 * [-1 1 1 -1];
    y = positions_mm(n,2) + element_y_mm/2 * [-1 -1 1 1];
    patch(x, y, corrected_delay_us(n), 'EdgeColor', 'w');
end

axis equal;
xlabel('x [mm]', 'FontSize', 14);
ylabel('y [mm]', 'FontSize', 14);
title('Aberration-corrected delays', 'FontSize', 17);
clim([0 delay_max]);
cb = colorbar;
cb.FontSize = 14;
ylabel(cb, 'Delay [\mus]', 'FontSize', 14);
grid on;

colormap(turbo);

if save_plots

    fig_name = sprintf('element_delays_steering_%gmm_%gmm_focus', x_steer_mm, z_focus_mm);

    savefig(fig_delays, fullfile([fig_name '.fig']));
    exportgraphics(fig_delays, fullfile([fig_name '.png']), 'Resolution', 300);

end