%% Trying to recreate transducer using source class
% Started from example_at_linear_array_transducer.m and then changed it

% Script to simulate 256-elements Frauenhofer transducer device, recording
% p_max on the whole grid and then plotting in 2 slices of interest to
% reproduce Figure 6 (except panel D) from device's paper. Therefore also 
% used to recreate 3.1, 3.7 and Supp. Fig. 6.10.
%
% NOTE: 
%
% author: Maylis Mueller
% date: 
% last update: 

clearvars;
addpath '/data/hu_muellerm/MATLAB Add-Ons/Collections/k-Wave/k-Wave'
addpath('/data/tu_muellerm/utils');

% =========================================================================
% USER-DEFINED LITERALS
% =========================================================================

% user input
save_plots = false;
add_delays = true;
sensor_axis = 3; % which plane to record on, 1 2 or 3 (x y or z)
get_FMHW = true;
FMHW = 'max_pressure'; % either "spe_freq" (uses extractAmpPhase) 
                        % or "max_pressure" (just uses p_max directly)

% computational parameters
t_end           = 160e-6;    % total compute time [s]
cfl             = 0.2;      % CFL number

%% medium parameters (for now water, would have to add skull and skin later)
c0              = 1480;     % sound speed [m/s] [water]
rho0            = 1000;     % density [kg/m^3] [water]

% assign medium properties
medium.sound_speed = c0;
medium.density = rho0;
medium.alpha_coeff = 2.17e-3;
medium.alpha_power = 2;

%% grid parameters (centered around 0 !!)
grid_size_x     = 80e-3;    % [m]
grid_size_y     = 80e-3;    % [m]
grid_size_z     = 150e-3;    % [m]

% calculate the grid spacing based on the PPW and F0
% dx = c0 / (ppw * source_f0);   % [m], 1.5e3/(3*5e5)= 1e-3 here
dx = 0.25e-3; % for now we just take a small dx, not based on calculation
dy = dx;
dz = dx;

%%
% =========================================================================
% CREATE THE GRID
% =========================================================================

% compute the size of the grid [grid points]
Nx = roundEven(grid_size_x / dx); % here 60
Ny = roundEven(grid_size_y / dx);
Nz = roundEven(grid_size_z / dx);

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% create the time array
kgrid.makeTime(c0, cfl, t_end); 
% here dt = cfl * dx / sound_speed 
% alternative:
% kgrid.makeTime(medium.sound_speed, cfl);
sample_freq = 1/kgrid.dt;


%%
% =========================================================================
% CREATE THE SOURCE => look at tus device GUI
% =========================================================================

%% source parameters
source_f0       = 4.8e5;      % source frequency [Hz]
% (initial values used purely out of theory
source_amp = 0.1005e6;
source_cycles   = 30;        % number of toneburst cycles

% treansducer geometry
row_counts = [6 10 12 14 16 16 18 18 18 18 18 18 16 16 14 12 10 6];   % custom aperture shape
element_width  = 3.9e-3;    % width [m]
element_height = element_width;    % Height of each element [m]
element_pitch  = 3.08e-3; % pitch [m]

% define focus position (in meters): cartesian distance from transducer center ;
% x and y in distance from transducer center
focal_point_user = [-20, -10, 40]*1e-3; 


%% ------------------------------------------------------------------------
% APPLY AXIAL FOCAL-DEPTH CALIBRATION
% -------------------------------------------------------------------------

% adding calculated delay
addpath('/data/tu_muellerm/delay_calibration/delay_calib_results');
load('delay_calibration.mat', 'focus_lengths', ...
    'corrected_sim_input');

if add_delays
    % focal_point_user is in metres, calibration data are in mm
    z_GUI_mm = focal_point_user(3) * 1e3;
    
    % Check that requested depth is within calibrated range
    if z_GUI_mm < min(focus_lengths) || z_GUI_mm > max(focus_lengths)
        error(['Requested focal depth %.1f mm is outside the calibrated ' ...
               'range [%.1f, %.1f] mm.'], ...
               z_GUI_mm, min(focus_lengths), max(focus_lengths));
    end
    
    % Convert GUI/requested z-position to corrected simulation input
    z_corrected_mm = interp1( ...
        focus_lengths, ...
        corrected_sim_input, ...
        z_GUI_mm, ...
        'pchip');
    z_corrected = z_corrected_mm * 1e-3;

else
    z_corrected = focal_point_user(3);
end

% x and y remain unchanged; only z is calibrated
focal_point_corrected = [ ...
    focal_point_user(1), ...
    focal_point_user(2), ...
    z_corrected];

% Convert to k-Wave grid coordinates
focal_point_grid = [ ...
    focal_point_corrected(1), ...
    focal_point_corrected(2), ...
    kgrid.z_vec(1) + focal_point_corrected(3)];


% repetition period between bursts
PRF_period = 15e-6;   % 15 microseconds

% number of bursts
num_bursts = 1;

[source, signal_info] = create_source(focal_point_grid, row_counts, ...
        element_height, element_width, element_pitch, ...
        source_amp, sample_freq, source_f0, c0, kgrid, ...
        'NumCycles', source_cycles);


%%
% --------------------
% SENSOR
% --------------------

% set sensor mask to record central plane
sensor.mask = ones(Nx, Ny, Nz);

% record the pressure
sensor.record = {'p_max'};

%%
% =========================================================================
% SIMULATION
% =========================================================================

input_args = {...
    %'PMLAlpha', PML_alpha, ...
   'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off'};

sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, ...
    input_args{:}, ...
    'DataCast', 'single', ...
    'PlotScale', [-1, 1] * source_amp, 'PlotSim', false);

p_max = reshape(sensor_data.p_max, [Nx, Ny, Nz]);

% --- axis-dependent grid vectors and labels ---

grid_vecs_2        = {kgrid.x_vec, kgrid.z_vec};
plane_dims_2       = [Nx, Nz];
axis_labels_2      = {'x-position [mm]', 'z-position [mm]'};
plane_label_2      = 'x–z plane';
prop_vec_shifted_2 = 1e3 * (kgrid.z_vec - kgrid.z_vec(1));  % along x axis
[~, focal_idx_2] = min(abs(kgrid.x_vec - focal_point_grid(2))); 
p_max_2 = squeeze(p_max(:, focal_idx_2, :));

grid_vecs_3        = {kgrid.x_vec, kgrid.y_vec};
plane_dims_3       = [Nx, Ny];
axis_labels_3      = {'x-position [mm]', 'y-position [mm]'};
plane_label_3      = 'x–y plane';
prop_vec_shifted_3 = 1e3 * kgrid.y_vec;                      % along y axis
[~, focal_idx_3] = min(abs(kgrid.z_vec - focal_point_grid(3))); 
p_max_3 = squeeze(p_max(:, :, focal_idx_3));

%%
% =========================================================================
% VISUALISATION
% =========================================================================

% x-z plane: plot along z at x=0
    axis_plot_label_2 = 'z-position [mm]';
    axis_name_2 = 'x';
% x-y plane: plot along y at x=0
    axis_plot_label_3 = 'x-position [mm]';
    axis_name_3 = 'y';
        
% -------------------------------------------------------------------------
% labels and filenames
% -------------------------------------------------------------------------

plot_info = sprintf( ...
    ['%.1fmm element size | %d bursts | %d cycles | ' ...
    'Focus = [%.1f %.1f %.1f] mm'], ...
    element_width*1e3, num_bursts, source_cycles, ...
    1e3 * focal_point_user(1), ...
    1e3 * focal_point_user(2), ...
    1e3 * focal_point_user(3));

filename_info = sprintf( ...
    '%.1fmm_el_size_%dbursts_%dcycles_focus_%0.0f_%0.0f_%0.0fmm_axis%d', ...
    element_width*1e3, num_bursts, source_cycles, ...
    1e3 * focal_point_user(1), ...
    1e3 * focal_point_user(2), ...
    1e3 * focal_point_user(3), ...
    sensor_axis);

output_folder = fullfile(pwd, 'simulation_plots/Figure_6');

% -------------------------------------------------------------------------
% plot pressure field in sensor plane
% -------------------------------------------------------------------------

% first in MPa
fig2 = figure;
imagesc( ...
    1e3 * grid_vecs_2{1}, ...
    prop_vec_shifted_2, ...
    (1e-6 * p_max_2).');
xlabel(axis_labels_2{1});
ylabel(axis_labels_2{2});
axis image;
title(sprintf('Pressure Field (%s)\n%s\n', plane_label_2, plot_info));
cb = colorbar;
title(cb, '[MPa]');


fig3 = figure;
imagesc( ...
    1e3 * grid_vecs_3{1}, ...
    prop_vec_shifted_3, ...
    (1e-6 * p_max_3).');
xlabel(axis_labels_3{1});
ylabel(axis_labels_3{2});
axis image;
title(sprintf('Pressure Field (%s)\n%s\n', plane_label_3, plot_info));
cb = colorbar;
title(cb, '[MPa]');

% then in dB

% paper-style colormap (Fig. 6): black -> blue -> green -> yellow -> red -> white
cmap_pos = [0.00 0.18 0.35 0.50 0.63 0.75 0.85 0.885 0.905 0.945 1.00];
cmap_rgb = [0     0     0   
            0     0     0.35
            0     0.05  0.75
            0.10  0.35  1.00
            0     0.35  0.55
            0.05  0.85  0.10
            1     1     0   
            1     0.3   0   
            1     0     0   
            1     0     0.05
            1     1     1   ];
cmap_db  = interp1(cmap_pos, cmap_rgb, linspace(0, 1, 256).', 'linear');
db_ticks = [-40 -20 -10 -6 -3 0];

db_floor = -40;

% axis 2
p_max_norm_2 = p_max_2 / max(p_max_2(:));   % normalize amplitude
p_max_db_2 = 20 * log10(p_max_norm_2);      % convert to dB

% avoid extremely low values dominating the color scale
p_max_db_2(p_max_db_2 < db_floor) = db_floor;

fig_db2 = figure;
imagesc( ...
    1e3 * grid_vecs_2{1}, ...
    prop_vec_shifted_2, ...
    p_max_db_2.');
xlabel(axis_labels_2{1});
ylabel(axis_labels_2{2});
title(sprintf('Pressure Field (dB)\n%s\n', plot_info));
axis image;
cb = colorbar;
title(cb, '[dB]');
cb.Ticks = db_ticks;

colormap(fig_db2, cmap_db);
clim([db_floor 0]);

% axis 3
p_max_norm_3 = p_max_3 / max(p_max_3(:));   % normalize amplitude
p_max_db_3 = 20 * log10(p_max_norm_3);      % convert to dB

p_max_db_3(p_max_db_3 < db_floor) = db_floor;

fig_db3 = figure;
imagesc( ...
    1e3 * grid_vecs_3{1}, ...
    prop_vec_shifted_3, ...
    p_max_db_3.');
xlabel(axis_labels_3{1});
ylabel(axis_labels_3{2});
title(sprintf('Pressure Field (dB)\n%s\n', plot_info));
axis image;
cb = colorbar;
title(cb, '[dB]');
cb.Ticks = db_ticks;

colormap(fig_db3, cmap_db);
clim([db_floor 0]);

% -------------------------------------------------------------------------
% save
% -------------------------------------------------------------------------

if save_plots

    % create folder
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    % save 2D pressure distribution in MPa
    savefig(fig2, fullfile(output_folder, ...
        sprintf('PressureField_axis_2_%s.fig', filename_info))); % axis 2
    saveas(fig2, fullfile(output_folder, ...
        sprintf('PressureField_axis_3_%s.png', filename_info)));

    savefig(fig3, fullfile(output_folder, ...
        sprintf('PressureField_axis_3%s.fig', filename_info))); % axis 3
    saveas(fig3, fullfile(output_folder, ...
        sprintf('PressureField_axis_3_%s.png', filename_info)));

    % save 2D pressure distribution in dB
    savefig(fig_db2, fullfile(output_folder, ...
        sprintf('PressureField_dB_axis_2_%s.fig', filename_info)));
    saveas(fig_db2, fullfile(output_folder, ...
        sprintf('PressureField_dB_axis_2_%s.png', filename_info)));

    savefig(fig_db3, fullfile(output_folder, ...
        sprintf('PressureField_dB_axis_3_%s.fig', filename_info)));
    saveas(fig_db3, fullfile(output_folder, ...
        sprintf('PressureField_dB_axis_3_%s.png', filename_info)));
end


% plot source output


% -------------------------------------------------------------------------
% plot FMHW
% -------------------------------------------------------------------------

    
if get_FMHW
    TargetedFocus = {};
    MeasuredFocus = {};
    
    FWHM_x_all = [];
    FWHM_z_all = [];
    
    PeakIndex_x_all = [];
    PeakIndex_z_all = [];
    
    FWHM_x_peak_all = [];
    FWHM_z_peak_all = [];
    
    % Find center indices
    [~, cx] = min(abs(kgrid.x_vec - focal_point_grid(1)));
    [~, cz] = min(abs(kgrid.z_vec - focal_point_grid(3)));
    
    
    % -------------------------
    % X PROFILE (at focal z)
    % -------------------------
    amp_x = p_max_2(:, cz);
    amp_x = amp_x / max(amp_x);
    
    halfMax_x = 0.5;
    
    left_x  = find(amp_x >= halfMax_x, 1, 'first');  
    right_x = find(amp_x >= halfMax_x, 1, 'last');
    
    disp('Left x :') ; disp(left_x);
    disp('Right x :') ; disp(right_x);
    
    fwhm_x = (right_x - left_x) * dx * 1e3; % mm
    
    % -------------------------
    % Z PROFILE (at focal x)
    % -------------------------
    amp_z = p_max_2(cx, :);
    amp_z = amp_z / max(amp_z);
    
    halfMax_z = 0.5;
    
    left_z  = find(amp_z >= halfMax_z, 1, 'first');
    right_z = find(amp_z >= halfMax_z, 1, 'last');
    
    disp('Left z :') ; disp(left_z);
    disp('Right z :') ; disp(right_z);
    
    fwhm_z = (right_z - left_z) * dz * 1e3; % mm
    
    % find peak location in the 2D field
    [~, maxIdx] = max(p_max_2(:));
    [ix_peak, iz_peak] = ind2sub(size(p_max_2), maxIdx);
    
    x_peak_mm = kgrid.x_vec(ix_peak) * 1e3;
    z_peak_mm = (kgrid.z_vec(iz_peak) + grid_size_z/2) * 1e3;

    % find peak location in OTHER 2D field
    [~, maxIdx] = max(p_max_3(:));
    [ix_peak, iy_peak] = ind2sub(size(p_max_3), maxIdx);
    
    y_peak_mm = kgrid.x_vec(iy_peak) * 1e3;
    
    % -------------------------
    % X PROFILE (at max z)
    % -------------------------
    amp_x_max_z = p_max_2(:, iz_peak);
    amp_x_max_z = amp_x_max_z / max(amp_x_max_z);
    
    halfMax_x_max_z = 0.5;
    
    left_x_max_z  = find(amp_x_max_z >= halfMax_x_max_z, 1, 'first');  
    right_x_max_z = find(amp_x_max_z >= halfMax_x_max_z, 1, 'last');
    
    disp('Left x at max z:') ; disp(left_x_max_z);
    disp('Right x at max z:') ; disp(right_x_max_z);
    
    fwhm_x_max_z = (right_x_max_z - left_x_max_z) * dx * 1e3; % mm
    
    % -------------------------
    % Z PROFILE (at max x)
    % -------------------------
    amp_z_max_x = p_max_2(ix_peak, :);
    amp_z_max_x = amp_z_max_x / max(amp_z_max_x);
    
    halfMax_z_max_x = 0.5;
    
    left_z_max_x  = find(amp_z_max_x >= halfMax_z_max_x, 1, 'first');
    right_z_max_x = find(amp_z_max_x >= halfMax_z_max_x, 1, 'last');
    
    disp('Left z at max x:') ; disp(left_z_max_x);
    disp('Right z at max x:') ; disp(right_z_max_x);
    
    fwhm_z_max_x = (right_z_max_x - left_z_max_x) * dz * 1e3; % mm
    
    
    % -------------------------
    % ADD VALUES TO LISTS
    % -------------------------
    
    % targeted focus string
    TargetedFocus{end+1} = sprintf( ...
        '%.0f/%.0f/%.0f', ...
        focal_point_user(1)*1e3, ...
        focal_point_user(2)*1e3, ...
        focal_point_user(3)*1e3);
    
    % measured focus string
    MeasuredFocus{end+1} = sprintf( ...
        '%.1f/%.1f', ...
        x_peak_mm, ...
        z_peak_mm);
    
    % numeric metrics
    FWHM_x_all(end+1) = fwhm_x;
    FWHM_z_all(end+1) = fwhm_z;
    
    FWHM_x_peak_all(end+1) = fwhm_x_max_z;
    FWHM_z_peak_all(end+1) = fwhm_z_max_x;
    
    % -------------------------
    % PLOT FWHM PROFILES
    % -------------------------
    
    fig_fmhw = figure;
    
    % X profile
    subplot(1,2,1)
    plot(kgrid.x_vec * 1e3, amp_x_max_z, 'LineWidth', 1.5); hold on;
    yline(0.5, '--r', 'Half Max');
    xline(kgrid.x_vec(left_x_max_z) * 1e3, '--k');
    xline(kgrid.x_vec(right_x_max_z) * 1e3, '--k');
    
    plot(kgrid.x_vec(left_x_max_z:right_x_max_z) * 1e3, ...
         amp_x_max_z(left_x_max_z:right_x_max_z), ...
         'LineWidth', 3);
    
    xlabel('x-position [mm]');
    ylabel('Normalized amplitude');
    title(sprintf('X Profile @ z = %.1f mm\nFWHM = %.2f mm', ...
        z_peak_mm, fwhm_x_max_z));
    grid on;
    ylim([0 1.1]);
    
    % Z profile
    subplot(1,2,2)
    plot(prop_vec_shifted_2, amp_z_max_x, 'LineWidth', 1.5); hold on;
    yline(0.5, '--r', 'Half Max');
    xline(prop_vec_shifted_2(left_z_max_x), '--k');
    xline(prop_vec_shifted_2(right_z_max_x), '--k');
    
    plot(prop_vec_shifted_2(left_z_max_x:right_z_max_x), ...
         amp_z_max_x(left_z_max_x:right_z_max_x), ...
         'LineWidth', 3);
    
    xlabel('z-position [mm]');
    ylabel('Normalized amplitude');
    title(sprintf('Z Profile @ x = %.1f mm\nFWHM = %.2f mm', ...
        x_peak_mm, fwhm_z_max_x));
    grid on;
    ylim([0 1.1]);
    
    sgtitle(sprintf( ...
        'FWHM Profiles\nInput Focus = [%.1f, %.1f, %.1f] mm ; Actual Focus = [%.1f, %.1f, %.1f] mm', ...
        focal_point_user(1)*1e3, ...
        focal_point_user(2)*1e3, ...
        focal_point_user(3)*1e3, ...
        x_peak_mm, ...
        y_peak_mm, ...
        z_peak_mm));
    
    Table1 = table( ...
        TargetedFocus(:), ...
        FWHM_x_all(:), ...
        FWHM_z_all(:), ...
        FWHM_x_peak_all(:), ...
        FWHM_z_peak_all(:), ...
        MeasuredFocus(:), ...
        'VariableNames', { ...
        'TargetedFocus', ...
        'FWHM_x_centerline_mm', ...
        'FWHM_z_centerline_mm', ...
        'FWHM_x_at_peak_z_mm', ...
        'FWHM_z_at_peak_x_mm', ...
        'MeasuredFocus_x_z_mm', ...
        });
    
    if save_plots
        savefig(fig_fmhw, fullfile(output_folder, ...
    sprintf('FMHW_%s_%s.fig', FMHW, filename_info)));
    
        saveas(fig_fmhw, fullfile(output_folder, ...
    sprintf('FMHW_%s_%s.png', FMHW, filename_info)));
    end
end


% save(sprintf('pressure_PML_alpha_%d_size_%d', PML_alpha, PML_size), 'pressure')