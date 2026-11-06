function [signal, info] = create_signal(sample_freq, source_f0, source_amp, Nt, varargin)
%CREATE_SIGNAL Create the temporal excitation signal for the TRUST source.
%
% The signal is constructed hierarchically as:
%
%   carrier cycles
%       -> pulse
%       -> pulse pause
%       -> pulse train
%       -> pulse-train pause
%       -> repeated pulse trains
%
% The pulse duration can be specified either by:
%
%   'NumCycles', N
%
% or
%
%   'PulseDuration', T
%
% but not both.
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS
% -------------------------------------------------------------------------
%
%   sample_freq    sampling frequency [Hz]
%   source_f0      carrier frequency [Hz]
%   source_amp     source pressure amplitude [Pa]
%   Nt             number of time samples to generate
%
% -------------------------------------------------------------------------
% NAME-VALUE INPUTS
% -------------------------------------------------------------------------
%
%   'NumCycles'             number of carrier cycles in one pulse
%
%   'PulseDuration'         requested pulse duration [s]
%
%   'PulsePause'            pause following each pulse [s]
%                           default = 0
%
%   'PRFPeriod'             alternatively specify pulse start-to-start
%                           period directly [s]
%
%   'NumPulsesPerTrain'     number of pulses in one pulse train
%                           default = 1
%
%   'PulseTrainPause'       pause following one pulse train [s]
%                           default = 0
%
%   'NumTrains'             number of repeated pulse trains
%                           default = 1
%
%   'RampUpDuration'        lower boundary of toneBurst enveloppe
%                           default =  Gaussian
%
%   'RampDownDuration'      upper boundary of toneBurst enveloppe
%                           default =  Gaussian
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
%
%   signal    1 x Nt temporal pressure signal [Pa]
%
%   info      structure containing the actual timing parameters used
%
% -------------------------------------------------------------------------
%   author: Maylis Mueller

% =========================================================================
% PARSE INPUTS
% =========================================================================

p = inputParser;

addParameter(p, 'NumCycles', []);
addParameter(p, 'PulseDuration', []);

addParameter(p, 'PulsePause', []);
addParameter(p, 'PRFPeriod', []);

addParameter(p, 'NumPulsesPerTrain', 1);
addParameter(p, 'PulseTrainPause', 0);
addParameter(p, 'NumTrains', 1);

% Envelope
addParameter(p, 'RampUpDuration', []);     % [ms]
addParameter(p, 'RampDownDuration', []);   % [ms]

parse(p, varargin{:});

num_cycles          = p.Results.NumCycles;
pulse_duration      = p.Results.PulseDuration;
pulse_pause         = p.Results.PulsePause;
PRF_period          = p.Results.PRFPeriod;
num_pulses_per_train = p.Results.NumPulsesPerTrain;
pulse_train_pause   = p.Results.PulseTrainPause;
num_trains          = p.Results.NumTrains;
ramp_up_ms   = p.Results.RampUpDuration;
% ramp_down_ms = p.Results.RampDownDuration;

ramp_down_ms = ramp_up_ms;

% =========================================================================
% CHECK PULSE-DURATION INPUT
% =========================================================================

has_cycles   = ~isempty(num_cycles);
has_duration = ~isempty(pulse_duration);

if ~has_cycles && ~has_duration
    error('Specify either ''NumCycles'' or ''PulseDuration''.');

elseif has_cycles && has_duration
    if abs(pulse_duration - num_cycles/source_f0) > 1e-12
        error(['If both ''NumCycles'' and ''PulseDuration'' are given, ' ...
               'they must correspond to the same pulse duration.']);
    end
end


% =========================================================================
% CONVERT PULSE DURATION <-> NUMBER OF CYCLES
% =========================================================================

if has_duration

    requested_cycles = pulse_duration * source_f0;

    num_cycles = round(requested_cycles);

    if abs(num_cycles - requested_cycles) > 1e-10

        warning(['Requested pulse duration corresponds to %.3f cycles. ' ...
                 'Using %d complete cycles instead.'], ...
                 requested_cycles, num_cycles);
    end

end

if num_cycles <= 0 || num_cycles ~= round(num_cycles)
    error('NumCycles must be a positive integer.');
end

% physical duration represented by the chosen number of cycles
pulse_duration = num_cycles / source_f0;


% =========================================================================
% DEFINE ENVELOPE
% =========================================================================

has_ramp_up   = ~isempty(ramp_up_ms);
has_ramp_down = ~isempty(ramp_down_ms); % even thougb in device ramp_up always equal ramp_down

% Default: Gaussian envelope
if ~has_ramp_up && ~has_ramp_down || (ramp_up_ms == 0 && ramp_down_ms == 0)

    envelope = 'Gaussian';

    ramp_up_cycles   = [];
    ramp_down_cycles = [];

% User-defined ramp-up and ramp-down
elseif has_ramp_up && has_ramp_down && ramp_up_ms ~= 0 && ramp_down_ms ~= 0

    % Check input
    if ramp_up_ms < 0 || ramp_down_ms < 0
        error('Ramp-up and ramp-down durations must be positive.');
    end

    % Convert ramp durations [ms] to number of carrier cycles
    ramp_up_cycles = ...
        ramp_up_ms * 1e-3 * source_f0;

    ramp_down_cycles = ...
        ramp_down_ms * 1e-3 * source_f0;

    % Check that ramps fit inside the pulse
    if ramp_up_cycles + ramp_down_cycles > num_cycles

        error(['Ramp durations are too long for the pulse.\n' ...
               'Pulse duration     = %.3f ms\n' ...
               'Ramp-up duration   = %.3f ms\n' ...
               'Ramp-down duration = %.3f ms'], ...
               num_cycles / source_f0 * 1e3, ...
               ramp_up_ms, ...
               ramp_down_ms);
    end

    % k-Wave expects ramp lengths in number of cycles
    envelope = [ramp_up_cycles, ramp_down_cycles];

else

    error(['Specify both ''RampUpDuration'' and ' ...
           '''RampDownDuration'', or neither.']);

end


% =========================================================================
% CREATE ONE ULTRASOUND PULSE
% =========================================================================

single_pulse = toneBurst( ...
    sample_freq, ...
    source_f0, ...
    num_cycles, ...
    'Envelope', envelope);

% single_pulse = single_pulse(:).';

% Apply acoustic source amplitude
single_pulse = source_amp .* single_pulse;

N_pulse = length(single_pulse);


% =========================================================================
% DEFINE PULSE REPETITION TIMING
% =========================================================================

has_pause = ~isempty(pulse_pause);
has_PRF   = ~isempty(PRF_period);

if has_pause && has_PRF
    N_pause_requested = round(pulse_pause * sample_freq);
    N_pause_PRF = round(PRF_period * sample_freq) - N_pulse;

    if N_pause_requested ~= N_pause_PRF
        error(['Specify either ''PulsePause'' or ''PRFPeriod'', ' ...
               'not both, OR have them match.']);
    end
end

if has_PRF

    % preserve the requested start-to-start repetition period
    N_pulse_period = round(PRF_period * sample_freq);

    N_pulse_pause = N_pulse_period - N_pulse;

    if N_pulse_pause < 0
        error(['PRF period is shorter than the pulse itself.\n' ...
               'Pulse duration = %.3f us\n' ...
               'PRF period     = %.3f us'], ...
               N_pulse / sample_freq * 1e6, ...
               PRF_period * 1e6);
    end

else

    if isempty(pulse_pause)
        pulse_pause = 0;
    end

    N_pulse_pause = round(pulse_pause * sample_freq);

    if N_pulse_pause < 0
        error('PulsePause must be >= 0.');
    end

    N_pulse_period = N_pulse + N_pulse_pause;

end


% =========================================================================
% DEFINE PULSE TRAIN TIMING
% =========================================================================

if num_pulses_per_train < 1 || ...
        num_pulses_per_train ~= round(num_pulses_per_train)

    error('NumPulsesPerTrain must be a positive integer.');
end

if num_trains < 1 || num_trains ~= round(num_trains)
    error('NumTrains must be a positive integer.');
end

N_train_pause = round(pulse_train_pause * sample_freq);

if N_train_pause < 0
    error('PulseTrainPause must be >= 0.');
end

% One pulse "iteration" is:
%
%      pulse + pulse pause
%
% This corresponds to the TRUST GUI terminology.
N_train = num_pulses_per_train * N_pulse_period;

% Complete train-to-train repetition period
N_train_period = N_train + N_train_pause;


% =========================================================================
% BUILD SIGNAL
% =========================================================================

signal = zeros(1, Nt);

for train_idx = 0:num_trains-1

    train_start = train_idx * N_train_period + 1;

    % Nothing else will fit in the requested simulation time
    if train_start > Nt
        break;
    end

    for pulse_idx = 0:num_pulses_per_train-1

        start_idx = ...
            train_start + pulse_idx * N_pulse_period;

        if start_idx > Nt
            break;
        end

        stop_idx = min( ...
            start_idx + N_pulse - 1, ...
            Nt);

        N_copy = stop_idx - start_idx + 1;

        signal(start_idx:stop_idx) = ...
            signal(start_idx:stop_idx) + ...
            single_pulse(1:N_copy);

    end
end


% =========================================================================
% OUTPUT INFORMATION
% =========================================================================

info.num_cycles = num_cycles;

info.pulse_duration = ...
    N_pulse / sample_freq;

info.pulse_pause = ...
    N_pulse_pause / sample_freq;

info.pulse_period = ...
    N_pulse_period / sample_freq;

info.PRF = ...
    1 / info.pulse_period;

info.num_pulses_per_train = ...
    num_pulses_per_train;

info.pulse_train_duration = ...
    N_train / sample_freq;

info.pulse_train_pause = ...
    N_train_pause / sample_freq;

info.pulse_train_period = ...
    N_train_period / sample_freq;

info.num_trains = num_trains;

info.total_sequence_duration = ...
    num_trains * N_train_period / sample_freq;

info.pulse_duty_cycle = ...
    100 * N_pulse / N_pulse_period;

info.total_duty_cycle = ...
    100 * ...
    (num_pulses_per_train * N_pulse) / ...
    N_train_period;

% enveloppe info
if ischar(envelope)

    info.envelope = envelope;
    info.ramp_up_ms = [];
    info.ramp_down_ms = [];

else

    info.envelope = 'Cosine taper';

    info.ramp_up_ms = ramp_up_ms;
    info.ramp_down_ms = ramp_down_ms;

    info.ramp_up_cycles = ramp_up_cycles;
    info.ramp_down_cycles = ramp_down_cycles;

end

end