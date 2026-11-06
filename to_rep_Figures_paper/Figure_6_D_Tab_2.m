%% Trying to recreate Figure 6.D from TUS device paper

% Script to simulate 256-elements Frauenhofer transducer device and Figure
% 6.D and Table 2 of device's paper. Therefore also used to recreate 3.1.D,
% 3.7.D and Supp. Fig. 6.10.D.
%
% NOTE: 
%
% author: Maylis Mueller
% date: 
% last update: 

clearvars;

addpath '/data/hu_muellerm/MATLAB Add-Ons/Collections/k-Wave/k-Wave'
addpath '/data/tu_muellerm/utils/'

addpath('/data/tu_muellerm/delay_calibration/delay_calib_results');
load('delay_calibration.mat', 'focus_lengths', ...
    'corrected_sim_input');


% =========================================================================
% USER-DEFINED LITERALS
% =========================================================================

% user input
save_plots = false;
add_delays = false;
sensor_axis = 3; % which plane to record on, 1 2 or 3 (x y or z)

% computational parameters
t_end           = 160e-6;    % total compute time [s]
cfl             = 0.3;      % CFL number

% medium parameters (for now water, would have to add skull and skin later)
c0              = 1480;     % sound speed [m/s] [water]
rho0            = 1000;     % density [kg/m^3] [water]

% assign medium properties
medium.sound_speed = c0;
medium.density = rho0;

% grid parameters (centered around 0 !!)
grid_size_x     = 350e-3;    % [m]
grid_size_y     = 80e-3;    % [m]
grid_size_z     = 150e-3;    % [m]

% calculate the grid spacing based on the PPW and F0
% dx = c0 / (ppw * source_f0);   % [m], 1.5e3/(3*5e5)= 1e-3 here
dx = 0.5e-3; % for now we just take a small dx, not based on calculation
dy = dx;
dz = dx;

% source parameters
source_f0       = 4.8e5;      % source frequency [Hz]
source_amp      = 1.5e6;      % source pressure [Pa] [1/2 Pa max]
source_cycles   = 30;        % number of toneburst cycles

% treansducer geometry
element_width  = 3.9e-3;    % width [m]
element_height = 3.9e-3;    % Height of each element [m]
element_pitch  = 3.08e-3; % pitch [m]

% transducer position
% translation     = [10e-3, 10e-3, 10e-3]; % [m]
% rotation        = [0, 0, 0]; % [deg]

% Row definition (disc-like shape, symmetric)
row_counts = [6 10 12 14 16 16 18 18 18 18 18 18 16 16 14 12 10 6];   % custom aperture shape

% define focus position (in meters)
focal_point_user = [0, 0, 100] * 1e-3; 
% cartesian distance from transducer center ; x and y in distance from
% transducer center


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
% here dt = cfl * dx / sound_speed = 3.3e-7 and total time is 35e-6 => 106
% time points
% alternative
% kgrid.makeTime(medium.sound_speed);
sample_freq = 1/kgrid.dt;


% =========================================================================
% STEERING ANGLE SWEEP FOR FIGURE D
% =========================================================================

steering_angles_deg = [0, 5.7, 11, 17, 22];  % degrees, matches figure 6
% steering_angles_deg = [0, 5.7, 11.3, 16.7, 21.8, 26.5];  % degrees, matches table 2

focal_depth_z       = 100e-3;                  % fixed depth [m]
colors = lines(length(steering_angles_deg));

num_elements = sum(row_counts);

time_delays = zeros(1, num_elements);

% ------------------------------------------------------------------------
% APPLY AXIAL FOCAL-DEPTH CALIBRATION
% -------------------------------------------------------------------------

% adding calculated delay
addpath('/data/tu_muellerm/delay_calibration');
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

input_args = { ...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off'};

% set sensor mask to record central plane
sensor.mask = zeros(Nx, Ny, Nz);

[~, focal_idx] = min(abs(kgrid.z_vec - focal_point_grid(3)));
    sensor.mask(:, :, focal_idx) = 1;

% record the pressure
sensor.record = {'p_max', 'p'};

fig_D = figure; hold on;

for s = 1:length(steering_angles_deg)

    theta = steering_angles_deg(s);

    focal_point_grid(1) = focal_point_user(3) * tand(theta);
    disp("Focal point grid :");
    disp(focal_point_grid);

    [source, signal_info] = create_source(focal_point_grid, row_counts, ...
        element_height, element_width, element_pitch, ...
        source_amp, sample_freq, source_f0, c0, kgrid, ...
        'NumCycles', source_cycles);


    % Run simulation
    sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, ...
        input_args{:}, 'DataCast', 'single', ...
        'PlotScale', [-1,1]*source_amp);

    % extract amplitude
    amp = extractAmpPhase(sensor_data.p, 1/kgrid.dt, source_f0, ...
        'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
    amp = reshape(amp, [Nx, Ny]);

    % Profile along x at y=0
    [~, cy] = min(abs(kgrid.y_vec));
    amp_x = amp(:, cy);
    amp_x_norm = amp_x / max(abs(amp_x));
    amp_x_dB   = 20*log10(amp_x_norm + eps); % avoids log(0)
    % amp_x_dB(amp_x_dB < -60) = -60;

    amp_x_dB(amp_x_dB < -100) = -100;
    ylim([-100 0])
    
    angles_deg = atan2d(kgrid.x_vec, focal_depth_z);

    plot(angles_deg, amp_x_dB, 'Color', colors(s,:), ...
        'DisplayName', sprintf('%.1f°', theta));
end
output_folder = fullfile(pwd, 'simulation_plots');

xlabel('Angle [°]');
ylabel('Pressure [dB]');
title('Beam Pattern at z = 100 mm');
xlim([-60 60]);
ylim([-60 0]);
legend show;
grid on;

% =========================================================================
% EXTRACT SECOND-HIGHEST LOCAL PEAK FOR EACH CURVE
% =========================================================================

ax = findobj(fig_D, 'Type', 'axes');
lines = findobj(ax, 'Type', 'line');
lines = flip(lines);   % restore original plotting order

n_curves = numel(lines);

second_peak_x  = nan(n_curves, 1);
second_peak_dB = nan(n_curves, 1);
curve_name     = strings(n_curves, 1);

for i = 1:n_curves

    x = lines(i).XData(:);
    y = lines(i).YData(:);

    curve_name(i) = string(lines(i).DisplayName);

    % Find local maxima, ordered from highest to lowest
    [pks, locs] = findpeaks(y, x, 'SortStr', 'descend');

    % Second-highest local maximum
    if i >=3
        if numel(pks) >= 2
            second_peak_dB(i) = pks(2);
            second_peak_x(i)  = locs(2);
        end
    elseif i == 1
        if numel(pks) >= 2
            second_peak_dB(i) = pks(9);
            second_peak_x(i)  = locs(9);
        end
    else
        if numel(pks) >= 2
            second_peak_dB(i) = pks(5);
            second_peak_x(i)  = locs(5);
        end
    end
end

% Display results
results = table(curve_name, second_peak_x, second_peak_dB, ...
    'VariableNames', {'SteeringAngle', 'PeakAngle_deg', 'Pressure_dB'});

disp(results);

% Save figure
savefig(fig_D, fullfile(output_folder, 'Figure_6_D.fig'));