%% Create k-Wave animation that displays pulsating "TUS"
% TUS_anim.m
%
% The first thing I did using k-Wave!
%
% author: Maylis Mueller

% Grid size
Nx = 200;
Ny = 200;
dx = 1e-3;
dy = 1e-3;

kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Create empty mask
source_mask = zeros(Nx, Ny);

% Shift values (centered)
shift_x = 55;
shift_y = -10;

% ----- Draw "T" -----
source_mask(20+shift_x:25+shift_x, 40+shift_y:80+shift_y) = 1;     
source_mask(20+shift_x:70+shift_x, 58+shift_y:62+shift_y) = 1;     

% ----- Draw "U" -----
source_mask(20+shift_x:70+shift_x, 90+shift_y:95+shift_y) = 1;     
source_mask(20+shift_x:70+shift_x, 120+shift_y:125+shift_y) = 1;   
source_mask(65+shift_x:70+shift_x, 90+shift_y:125+shift_y) = 1;    

% ----- Draw "S" -----
source_mask(20+shift_x:25+shift_x, 140+shift_y:180+shift_y) = 1;   
source_mask(20+shift_x:45+shift_x, 140+shift_y:145+shift_y) = 1;   
source_mask(40+shift_x:45+shift_x, 140+shift_y:180+shift_y) = 1;   
source_mask(45+shift_x:70+shift_x, 175+shift_y:180+shift_y) = 1;   
source_mask(65+shift_x:70+shift_x, 140+shift_y:180+shift_y) = 1;       

% Assign to source
source.p_mask = source_mask;

% Time array
kgrid.makeTime(1500);

% Simple sinusoidal source
freq = 0.5e6; % 0.5 MHz
% source.p = sin(2*pi*freq*kgrid.t_array);
source.p = toneBurst(1/kgrid.dt, freq, 5); % 5 cycles

% Medium properties
medium.sound_speed = 1500;

% Sensor (just to run simulation)
sensor.mask = ones(Nx, Ny);

% Run simulation
input_args = {'PMLInside', false, ...
              'PlotSim', true, ...
              'PlotFreq', 5, ...
              'DisplayMask', source.p_mask, ...
              'RecordMovie', false, ...
              'MovieName', 'TUS_wave4', ...
              'MovieType', 'mp4'};

sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});