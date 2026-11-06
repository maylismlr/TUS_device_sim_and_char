function pac = estimate_pac_phase_reversal(kgrid, medium_water, medium_skull, ...
    target_grid, row_counts, element_height, element_width, element_pitch, ...
    source_f0, source_amp, input_args, varargin)
%ESTIMATE_PAC_PHASE_REVERSAL Reciprocal two-run PAC estimation for TRUST.
%
% A virtual point source is placed at the intended target and the complete
% TRUST kWaveArray is used as a receiver. Acoustic reciprocity means that
% the target-to-element transfer function equals the element-to-target
% transfer function in the linear reciprocal model used here. Consequently,
% all element responses are obtained simultaneously.
%
% The simulation is run once in water and once with the skull:
%
%   delta_phase(n) = angle(H_skull(n) * conj(H_water(n)))
%   correction_phase(n) = -delta_phase(n)
%
% The water subtraction isolates the skull-induced phase so that the PAC
% delay can be added to the existing calibrated geometric focusing delays.
%
% This reciprocal virtual-source approach follows the numerical
% time-reversal principle described by Jing, Meral, and Clement
% (Phys. Med. Biol. 57, 901-917, 2012). Here, the skull-induced
% phase shift at the carrier frequency is extracted relative to water
% and converted to an element-wise phase-delay correction.
%
% target_grid uses the same centred physical coordinates as kgrid [m].
% medium_water and medium_skull must be defined on the same kgrid.
%
% Name-value options:
%   'NumCycles'       calibration burst length (default 30)
%   'SaveFile'        output MAT file
%                     (default 'pac_TRUST_reciprocal.mat')
%   'PlotProgress'       plot the resulting correction (default true)
%   'SaveProgressPlot'   save the progress plot in the current directory
%                        (default false)
%   'SaveWaveforms'      store the 256 received time traces (default false)
%
%   author: Maylis Mueller

p = inputParser;
addParameter(p, 'NumCycles', 30);
addParameter(p, 'SaveFile', 'pac_TRUST_reciprocal.mat');
addParameter(p, 'PlotProgress', true);
addParameter(p, 'SaveProgressPlot', false);
addParameter(p, 'SaveWaveforms', false);
parse(p, varargin{:});
opt = p.Results;
z_pos = kgrid.z_vec(1);
[karray, positions] = build_karray(row_counts, element_height, ...
    element_width, element_pitch, z_pos);
num_elements = numel(karray.elements);

% Snap the requested target to the nearest computational grid point.
[~, ix] = min(abs(kgrid.x_vec - target_grid(1)));
[~, iy] = min(abs(kgrid.y_vec - target_grid(2)));
[~, iz] = min(abs(kgrid.z_vec - target_grid(3)));
target_index = [ix, iy, iz];
target_grid_actual = [kgrid.x_vec(ix), kgrid.y_vec(iy), kgrid.z_vec(iz)];

target_error = norm(target_grid_actual - target_grid);
if target_error > 0.51 * max([kgrid.dx, kgrid.dy, kgrid.dz])
    warning('PAC target was shifted by %.3f mm when mapped to the grid.', ...
        target_error * 1e3);
end

% One virtual pressure source at the intended target.
source.p_mask = false(kgrid.Nx, kgrid.Ny, kgrid.Nz);
source.p_mask(ix, iy, iz) = true;

sample_freq = 1 / kgrid.dt;
[source.p, signal_info] = create_signal(sample_freq, source_f0, ...
    source_amp, kgrid.Nt, 'NumCycles', opt.NumCycles);

% Use the same finite, off-grid TRUST elements as receivers. combineSensorData
% applies the appropriate BLI weights and returns one trace per element in
% the same ordering used by create_source_pac.
sensor.mask = karray.getArrayBinaryMask(kgrid);
sensor.record = {'p'};

source_component = component_at_frequency(source.p, kgrid.t_array, source_f0);
if abs(source_component) < eps
    error('The PAC source has no usable component at %.3f kHz.', ...
        source_f0 * 1e-3);
end

fprintf('PAC reciprocal simulation 1/2: target source, water medium\n');
data_water = kspaceFirstOrder3DG(kgrid, medium_water, source, sensor, input_args{:});
received_water = combine_element_data(karray, kgrid, data_water.p, ...
    num_elements);
H_water = component_at_frequency(received_water, kgrid.t_array, ...
    source_f0) ./ source_component;
clear data_water;

fprintf('PAC reciprocal simulation 2/2: target source, skull medium\n');

data_skull = kspaceFirstOrder3DG(kgrid, medium_skull, source, sensor, input_args{:});

received_skull = combine_element_data(karray, kgrid, data_skull.p, ...
    num_elements);
H_skull = component_at_frequency(received_skull, kgrid.t_array, ...
    source_f0) ./ source_component;
clear data_skull;

if any(~isfinite(H_water)) || any(~isfinite(H_skull))
    error('Non-finite PAC transfer-function values were obtained.');
end
if any(abs(H_water) < 100 * eps(max(abs(H_water))))
    warning(['At least one water response is very small; inspect the element ' ...
        'responses before using its phase correction.']);
end

delta_phase = angle(H_skull .* conj(H_water));
correction_phase = -delta_phase;

% A positive signal delay contributes -omega*tau to its carrier phase.
% Therefore tau = delta_phase/omega implements -delta_phase, modulo one
% carrier period. Removing a common delay does not change focusing.
correction_delay = mod(delta_phase, 2*pi) / (2*pi*source_f0);
correction_delay = correction_delay - min(correction_delay);

pac.method = 'reciprocal target source: skull minus water';
pac.delta_phase = delta_phase;
pac.correction_phase = correction_phase;
pac.correction_delay = correction_delay;
pac.H_water = H_water;
pac.H_skull = H_skull;
pac.transmission = abs(H_skull) ./ abs(H_water);
pac.positions = positions;
pac.target_user_requested = target_grid;
pac.target_user_requested(3) = target_grid(3)-kgrid.z_vec(1);
pac.target_grid = target_grid_actual;
pac.target_index = target_index;
pac.source_f0 = source_f0;
pac.num_forward_simulations_replaced = 2 * num_elements;
pac.num_reciprocal_simulations = 2;

if opt.SaveWaveforms
    save(opt.SaveFile, 'pac', 'received_water', ...
        'received_skull', 'positions', 'target_grid', 'target_grid_actual', ...
        'source_f0', 'source_amp', 'signal_info', '-v7.3');
else
    save(opt.SaveFile, 'pac', 'positions', ...
        'target_grid', 'target_grid_actual', 'source_f0', 'source_amp', ...
        'signal_info', '-v7.3');
end

if opt.PlotProgress
    fig_pac = figure;
    tiledlayout(3, 1);

    nexttile;
    plot(1:num_elements, rad2deg(delta_phase), '.');
    ylabel('\Delta\phi_{skull-water} [deg]');
    xlim([0 num_elements + 1]);
    grid on; box on;

    nexttile;
    plot(1:num_elements, correction_delay * 1e6, '.');
    ylabel('PAC delay [\mus]');
    xlim([0 num_elements + 1]);
    grid on; box on;

    nexttile;
    plot(1:num_elements, pac.transmission, '.');
    xlabel('TRUST element');
    ylabel('|H_{skull}|/|H_{water}|');
    xlim([0 num_elements + 1]);
    grid on; box on;

    if opt.SaveProgressPlot
        exportgraphics(fig_pac, fullfile(pwd, sprintf('PAC_progress_%.0f_%.0fmm.png', ...
            pac.target_user_requested(1)*1e3, pac.target_user_requested(3)*1e3)), 'Resolution', 300);
    
        savefig(fig_pac, fullfile(pwd, sprintf('PAC_progress_%.0f_%.0fmm.fig', ...
            pac.target_user_requested(1)*1e3, pac.target_user_requested(3)*1e3)));
    end
end
end

function element_data = combine_element_data(karray, kgrid, grid_data, ...
    num_elements)
    
    % Convert BLI grid-sensor traces to one receive trace per physical element.
    element_data = karray.combineSensorData(kgrid, grid_data);
    
    if size(element_data, 1) == num_elements
        return;
    elseif size(element_data, 2) == num_elements
        element_data = element_data.';
    else
        error('Could not identify the element dimension in combined sensor data.');
    end
end

function H = component_at_frequency(waveforms, t, f0)
% Complex Fourier coefficient at f0 using the absolute k-Wave time axis.
    x = double(waveforms);
    t = double(t(:));
    
    if isvector(x)
        x = x(:).';
    end
    if size(x, 2) ~= numel(t)
        if size(x, 1) == numel(t)
            x = x.';
        else
            error('Waveform length does not match kgrid.t_array.');
        end
    end
    
    x = x - mean(x, 2);
    carrier = exp(-1i * 2*pi*f0 .* t);
    H = x * carrier;
    H = H(:).';
end