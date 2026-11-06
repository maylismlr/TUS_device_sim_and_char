%% Trying to recreate transducer using source class
% Script to run the simulations with a muCT in the medium, then plot and 
% save pmax pressure field
%
% authors: Thomas Kirchner & Maylis Mueller

clearvars;
addpath '/data/tu_muellerm/utils/'
addpath '/data/hu_muellerm/MATLAB Add-Ons/Collections/k-Wave/k-Wave'

% =========================================================================
% USER-DEFINED LITERALS
% =========================================================================

path_medium = '/data/tu_muellerm_cloud/skulls_FUN/skulls_params/final/muCT.nrrd_0.25mm_params.mat';
tmp = load(path_medium);
seg = tmp.seg();

% computational parameters
t_end           = 120e-6;    % total compute time [s]
cfl             = 0.2;      % CFL number

% medium parameters
c0              = 1480;     % sound speed [m/s] [water]
rho0            = 1000;     % density [kg/m^3] [water]

% assign medium properties
medium.sound_speed = c0;
medium.density = rho0;

% medium.alpha_coeff = 1.10e-3;
% medium.alpha_power = 2.8;


% =========================================================================
% CREATE THE GRID
% =========================================================================

% grid parameters (centered around 0 !!)
grid_size_x     = 56e-3;    % [m]
grid_size_y     = 56e-3;    % [m]
grid_size_z     = 70e-3;    % [m]

% calculate the grid spacing based on the PPW and F0
dx = 0.25e-3; % for now we just take a small dx, not based on calculation
dy = dx;
dz = dx;

% compute the size of the grid [grid points]
Nx = roundEven(grid_size_x / dx);
Ny = roundEven(grid_size_y / dx);
Nz = roundEven(grid_size_z / dx);

medium.sound_speed = ones(Nx,Ny,Nz);
seg = seg(1:Nx, 1:Ny, 1:Nz);
medium.sound_speed(seg==0) = 1480.0;
medium.sound_speed(seg==1) = 3100.0;
medium.density = ones(Nx,Ny,Nz);
medium.density(seg==0) = 1000.0;
medium.density(seg==1) = 2000.0;
medium.alpha_coeff = ones(Nx,Ny,Nz);
medium.alpha_coeff(seg==0) = 2.17e-3;      % dB/(MHz^2 cm)
medium.alpha_coeff(seg==1) = 2.7;
medium.alpha_power = 2;

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% create the time array
kgrid.makeTime(3100.0, cfl, t_end); 
sample_freq = 1/kgrid.dt;

% define focus position (in m): cartesian distance from transducer center ;
% x and y in distance from transducer center
focal_point_user = [0, 0, 45]*1e-3; 

% now convert it in terms of the grid, so taking into account the fact that
% grid is centered around 0
focal_point_grid = [focal_point_user(1), ...
                    focal_point_user(2), ...
                    kgrid.z_vec(1) + focal_point_user(3)]; 


% =========================================================================
% CREATE THE SOURCE 
% =========================================================================

% source parameters
source_f0       = 4.8e5;     % source frequency [Hz]
source_volt = 30;
% source_amp      = 1.5e6;     % source pressure [Pa] [1/2 Pa max]
p_ref = trust_reference_data();
[source_amp, info] = trust_p0(source_volt, focal_point_user(3)*1e3, p_ref);
source_cycles   = 30;        % number of toneburst cycles

% treansducer geometry
element_width  = 3.9e-3;    % width [m]
element_height = element_width;    % Height of each element [m]
element_pitch  = 3.08e-3; % pitch [m]

% Row definition (disc-like shape, symmetric)
row_counts = [6 10 12 14 16 16 18 18 18 18 18 18 16 16 14 12 10 6];   % custom aperture shape


% repetition period between bursts
PRF_period = 15e-6;   % 15 microseconds

% number of bursts
num_bursts = 1;

[source, signal_info] = create_source(focal_point_grid, row_counts, ...
        element_height, element_width, element_pitch, ...
        source_amp, sample_freq, source_f0, c0, kgrid, ...
        'NumCycles', source_cycles);


% =========================================================================
% CREATE THE SENSOR 
% =========================================================================

% set sensor mask to record central plane
sensor.mask = ones(Nx, Ny, Nz);


% record max pressure only
sensor.record = {'p_max'};

%
% =========================================================================
% SIMULATION
% =========================================================================

input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off'};

sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, ...
    input_args{:}, ...
    'DataCast', 'single', ...
    'PlotScale', [-1, 1] * source_amp);

grid_vecs        = {kgrid.x_vec, kgrid.z_vec};
plane_dims       = [Nx, Nz];
axis_labels      = {'x-position [mm]', 'z-position [mm]'};
plane_label      = 'x–z plane';
prop_vec_shifted = 1e3 * (kgrid.z_vec-kgrid.z_vec(1));


% reshape p_max into the sensor plane
p_max = reshape(sensor_data.p_max, [Nx,Ny,Nz]);

save('/data/tu_muellerm_cloud/skulls_FUN/sensor_data/muCT_0.25_power2_p_max_0mm_45mm.mat', ...
 'p_max', '-v7.3')


% =========================================================================
% VISUALISATION
% =========================================================================

% plot pressure field in MPa in sensor plane
fig2 = figure;
imagesc( ...
    1e3 * grid_vecs{1}, ...
    prop_vec_shifted, ...
    (1e-6 * squeeze(p_max(:,112,:))).');
xlabel(axis_labels{1});
ylabel(axis_labels{2});
axis image;
title(sprintf('Pressure Field (%s)\n%s\n', plane_label));
cb = colorbar;
title(cb, '[MPa]');
