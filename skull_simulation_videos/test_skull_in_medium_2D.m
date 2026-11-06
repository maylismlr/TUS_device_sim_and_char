% Trying out adding the skull in medium in my own way !

clearvars;
addpath '/data/hu_muellerm/MATLAB Add-Ons/Collections/k-Wave/k-Wave'
addpath 'utils'

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% computational parameters
ppw             = 10;        % number of points per wavelength
t_end           = 30e-6;    % total compute time [s]
cfl             = 0.3;      % CFL number


%% source parameters
source_f0       = 4.8e5;      % source frequency [Hz]
source_amp      = 1.5e6;      % source pressure [Pa] [1/2 Pa max]
source_cycles   = 30;        % number of toneburst cycles


%% grid parameters (centered around 0 !!)
grid_size_x     = 50e-3;    % [m]
grid_size_y     = 50e-3;    % [m]

% calculate the grid spacing based on the PPW and F0
% dx = c0 / (ppw * source_f0);   % [m], 1.5e3/(3*5e5)= 1e-3 here
dx = 0.1e-3; % for now we just take a small dx, not based on calculation
dy = dx;


%%
% =========================================================================
% CREATE THE GRID
% =========================================================================

% compute the size of the grid [grid points]
Nx = roundEven(grid_size_x / dx); % here 60
Ny = roundEven(grid_size_y / dx);

kgrid = kWaveGrid(Nx, dx, Ny, dy);

% create the time array
c0 = 3100;
kgrid.makeTime(c0, cfl, t_end);

dt = kgrid.dt;
% here dt = cfl * dx / sound_speed = 3.3e-7 and total time is 35e-6 => 106
% time points
% alternative
% kgrid.makeTime(medium.sound_speed);
sample_freq = 1/kgrid.dt;
source_diam     = 20e-3;    % piston diameter [m]


% calculate the actual CFL and PPW
disp(['PPW = ' num2str(c0 / (dx * source_f0))]);
disp(['CFL = ' num2str(c0 * dt / dx)]);

% --------------------
% SOURCE
% --------------------

% create time varying source
tone_burst_cycles = 4;
source_sig = source_amp * toneBurst(1/kgrid.dt, source_f0, tone_burst_cycles);

% assign binary mask
source.p_mask = zeros(Nx, Ny);
source.p_mask(1, round(Ny/2 - source_diam/(2*dx)) : round(Ny/2 + source_diam/(2*dx))) = 1;

% assign source signals
source.p = source_sig;


% --------------------
% MEDIUM
% --------------------

path = '/data/tu_muellerm_cloud/skulls/halle_skull_FB_rotated.mat';
tmp = load(path);
slice_index = 256;
skull2D = squeeze(tmp.seg(:, slice_index, :));  % removes singleton dimension
skull2D = rot90(skull2D, 1);

[Nx_full, Ny_full] = size(skull2D);

x_start = round((Nx_full - Nx)/2);
y_start = round((Ny_full - Ny)/2);


skull2D = imresize(skull2D, [Nx Ny], 'nearest');
seg = skull2D;
% seg = skull2D(x_start:x_start+Nx-1, y_start:y_start+Ny-1);
          
% define the properties of the propagation medium    
medium.sound_speed = ones(Nx,Ny);
medium.sound_speed(seg==0) = 1480.0;
medium.sound_speed(seg==1) = 3100.0;
medium.density = ones(Nx,Ny);
medium.density(seg==0) = 1000.0;
medium.density(seg==1) = 2190.0;
medium.alpha_coeff = ones(Nx,Ny);
medium.alpha_coeff(seg==0) = 0.0;
medium.alpha_coeff(seg==1) = 2.7;
medium.alpha_power = 2;

kgrid.makeTime(medium.sound_speed, cfl, t_end); 


% --------------------
% SENSOR
% --------------------

% set sensor mask to record central plane
sensor.mask = zeros(Nx, Ny);
sensor.mask(:, :) = 1;

% record the pressure
sensor.record = {'p'};


% --------------------
% SIMULATION
% --------------------

display_mask = double(seg) + double(source.p_mask);

% set input options
input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', display_mask};
        
% MATLAB CPU
sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, ...
    input_args{:}, ...
    'DataCast', 'single', ...
    'PlotScale', [-1, 1]*source_amp, ...
    'RecordMovie', false, ...
    'MovieName', 'skull_in_medium_test_2', ...
    'MovieType', 'mp4');
        

% extract amplitude from the sensor data
amp = extractAmpPhase(sensor_data.p, 1/kgrid.dt, source_f0, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% reshape data
amp = reshape(amp, Nx, Ny);

% extract pressure on axis
amp_on_axis = amp(:, 1);

% define axis vectors for plotting
x_vec = kgrid.x_vec(:, :) - kgrid.x_vec(1);
y_vec = kgrid.y_vec - kgrid.y_vec(1);


% =========================================================================
% VISUALISATION
% =========================================================================

% plot the source mask (pml is outside the grid in this example)
figure;
imagesc(1e3 * (kgrid.y_vec - kgrid.y_vec(1)), 1e3 * kgrid.x_vec, source.p_mask);
xlabel('y [mm]');
ylabel('x [mm]');
axis image;
title('Source Mask');

% plot the pressure field 
figure;
imagesc(1e3 * y_vec, 1e3 * x_vec, amp);
xlabel('Lateral Position [mm]');
ylabel('Axial Position [mm]');
axis image;
title('Pressure Field');

% =========================================================================
% TIME SERIES
% =========================================================================

% NO ABS VALUE (neg pressure gets substracted from skull set to max
% pressure)
% pressure = sensor_data.p;
% pressure = reshape(pressure, 500, 500, length(kgrid.t_array));
% maxPressure = max(abs(pressure(:)));
% pressure_and_skull = double(seg)*maxPressure + double(pressure);

pressure = sensor_data.p;
pressure = reshape(pressure, 500, 500, length(kgrid.t_array));
maxPressure = max(abs(pressure(:)));

% Having use abs of values so that it truly is just putting the skull on
% top
mask = logical(seg);
pressure_and_skull = double(pressure);
pressure_and_skull = pressure_and_skull .* (1 - mask) + (abs(pressure_and_skull) + maxPressure) .* mask;
img.pixelData = pressure_and_skull;

% optional metadata
img.metaData.encoding = 'raw';

% nrrdwrite('2D_skull_pressure.nrrd', img);

sliceIndex = 1;

f = figure;
ax = axes(f);

% Display first slice
hImg = imagesc(pressure_and_skull(:,:,sliceIndex), 'Parent', ax);

axis image
colorbar

% plot_layout_cmap = flipud(bone(256));
% colormap(ax, plot_layout_cmap);
% colormap(ax, turbo)    

COLOR_MAP = getColorMap();

colormap(ax, COLOR_MAP);

% Increase contrast
% clim([0 max(pressure_and_skull(:))])  

clim([-maxPressure maxPressure])

sld = uicontrol('Style','slider', ...
    'Min',1,'Max',size(pressure_and_skull,3), ...
    'Value',sliceIndex, ...
    'SliderStep',[1/(size(pressure_and_skull,3)-1) 10/(size(pressure_and_skull,3)-1)], ...
    'Position',[150 20 300 20], ...
    'Callback', @(src,evt) updateSlice(src, hImg, pressure_and_skull));
