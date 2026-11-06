%% Temporal signal of the simulation paradigm
%
% Recreates Figure 7 of the device paper (so Figure 2.5 of the Master's
% thesis), aka system signal pattern.
% Reproduces GUI input parameters.
%
% author: Maylis Mueller
% date:
% last update:

clear; clc;
addpath '/data/hu_muellerm/MATLAB Add-Ons/Collections/k-Wave/k-Wave'
addpath '/data/tu_muellerm/utils/'


% -------------------------------------------------------------------------
% PARAMETERS
% -------------------------------------------------------------------------

source_f0   = 480e3;       % [Hz]
sample_freq = 50e6;         % [Hz] only for plotting

% Specify EITHER pulse duration OR number of cycles
pulse_duration = [];     % [s]
num_cycles     = 30;       % e.g. 2400; leave [] if using pulse_duration

pulse_pause          = 0.01e-3;     % [s]
num_pulses_per_train = 2;

pulse_train_pause = 0.05e-3;      % [s]
num_trains        = 2;

ramp_up_ms = 0.02;
ramp_down_ms = 0.04;

% Plot appearance
plot_amplitude = 0.3;       % amplitude of normalized signal [a.u.]

y_min = -0.4;               % y-axis minimum
y_max =  0.4;               % y-axis maximum

time_before = 50e-6;          % empty time before t = 0 [s]
extra_time  = 0.05;          % fraction of one train period shown after


% -------------------------------------------------------------------------
% PULSE DURATION / NUMBER OF CYCLES
% -------------------------------------------------------------------------

if ~isempty(pulse_duration) && ~isempty(num_cycles)
    error('Specify either pulse_duration OR num_cycles, not both.');
end

if isempty(pulse_duration) && isempty(num_cycles)
    error('Specify either pulse_duration OR num_cycles.');
end

if isempty(num_cycles)

    num_cycles = round(pulse_duration * source_f0);

elseif isempty(pulse_duration)

    pulse_duration = num_cycles / source_f0;

end


% -------------------------------------------------------------------------
% TIMING
% -------------------------------------------------------------------------

pulse_period = pulse_duration + pulse_pause;

PRF = 1 / pulse_period;

pulse_train_duration = ...
    num_pulses_per_train * pulse_period;

pulse_train_period = ...
    pulse_train_duration + pulse_train_pause;

sequence_duration = ...
    num_trains * pulse_train_period;

% Add some empty space after the last train for plotting
T_plot = sequence_duration + extra_time * pulse_train_period;

Nt = round(T_plot * sample_freq) + 1;


% -------------------------------------------------------------------------
% CREATE SIGNAL
% -------------------------------------------------------------------------

[signal, info] = create_signal( ...
    sample_freq, ...
    source_f0, ...
    1, ...                          % arbitrary amplitude
    Nt, ...
    'NumCycles', num_cycles, ...
    'PulsePause', pulse_pause, ...
    'NumPulsesPerTrain', num_pulses_per_train, ...
    'PulseTrainPause', pulse_train_pause, ...
    'NumTrains', num_trains, ...
    'RampUpDuration', ramp_up_ms, 'RampDownDuration', ramp_down_ms);


% -------------------------------------------------------------------------
% NORMALIZE FOR PLOTTING
% -------------------------------------------------------------------------

signal_plot = plot_amplitude .* signal ./ max(abs(signal));

t = (0:Nt-1) / sample_freq;


% -------------------------------------------------------------------------
% AUTOMATIC TIME UNIT
% -------------------------------------------------------------------------

if T_plot >= 0.1

    time_scale = 1;
    time_unit  = 's';

elseif T_plot >= 1e-3

    time_scale = 1e3;
    time_unit  = 'ms';

else

    time_scale = 1e6;
    time_unit  = '\mus';

end
time_scale = 1e6;
time_unit  = '\mus';

% -------------------------------------------------------------------------
% PLOT
% -------------------------------------------------------------------------

figure('GraphicsSmoothing', 'off');

plot( ...
    t * time_scale, ...
    signal_plot, ...
    'LineWidth', 0.5);

xlabel(sprintf('Time [%s]', time_unit));
ylabel('Amplitude [a.u.]');

title(sprintf(['Pulse duration %.3g %s | ' ...
               '%d cycles | %.f %s pulse pause \n %.f %s ramp up | %.f %s ramp down' ...
               '\n %d pulses/train | %.f %s pulse train pause'], ...
               pulse_duration * 1e6, time_unit,...
               num_cycles, ...
               pulse_pause*1e6, time_unit, ...
               ramp_up_ms*1e3, time_unit, ramp_down_ms*1e3, time_unit, ...
               num_pulses_per_train, ...
               pulse_train_pause*1e6, time_unit));

xlim([ ...
    -time_before * time_scale, ...
     T_plot * time_scale]);
%xticks([0 2000 4000])
ylim([y_min y_max]);

grid on;
box on;


% -------------------------------------------------------------------------
% DISPLAY PARAMETERS
% -------------------------------------------------------------------------

fprintf('\n');
fprintf('Carrier frequency     : %.1f kHz\n', source_f0/1e3);
fprintf('Number of cycles      : %d\n', num_cycles);
fprintf('Pulse duration        : %.3f ms\n', pulse_duration*1e3);
fprintf('Pulse pause           : %.3f ms\n', pulse_pause*1e3);
fprintf('Pulse period          : %.3f ms\n', pulse_period*1e3);
fprintf('PRF                   : %.2f Hz\n', PRF);
fprintf('Pulses per train      : %d\n', num_pulses_per_train);
fprintf('Pulse train duration  : %.3f ms\n', pulse_train_duration*1e3);
fprintf('Pulse train pause     : %.3f ms\n', pulse_train_pause*1e3);
fprintf('Pulse train period    : %.3f ms\n', pulse_train_period*1e3);
fprintf('Number of trains      : %d\n', num_trains);


% -------------------------------------------------------------------------
% SAVE
% -------------------------------------------------------------------------

save_dir = './simulation_plots/figure_7';

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% Descriptive filename based on signal parameters
fig_name = sprintf( ...
    ['Fig7_signal_f0_%gkHz_pulse_%gms_pause_%gms_' ...
     'nPulse_%d_trainPause_%gms_nTrain_%d_'], ...
    source_f0 / 1e3, ...
    pulse_duration * 1e3, ...
    pulse_pause * 1e3, ...
    num_pulses_per_train, ...
    pulse_train_pause * 1e3, ...
    num_trains);

% Replace decimal points to keep filenames clean
fig_name = strrep(fig_name, '.', 'p');

% Current figure
fig = gcf;

% Rendering settings
set(fig, 'Renderer', 'painters');
drawnow;

% Save editable MATLAB figure
savefig(fig, ...
    fullfile(save_dir, [fig_name '.fig']));

% Save high-resolution PNG
saveas(fig, ...
    fullfile(save_dir, [fig_name '.png']));