%% Report and plot input pressure voltage relationships
% plot_input_pressure.m
%
% Reports the amplitude calibration and the checks that support it.
%
%   1. p_eff is proportional to V at each measured depth
%   2. the normalised voltage response is the same at both depths
%   3. k(z), the source amplitude required at the reference drive
%   4. measured and simulated depth profiles, normalised at 50 mm
%   5. prediction of the two voltage sweeps from the depth sweep alone
%
% Plots Figures 2.7 and 3.5 from Master's thesis.
%
% Checks 1, 2 and 5 concern the reference data and the factorisation; none of
% them tests the acoustic model. Only 3 and 4 do.
%
%
% Author: Maylis Mueller

clear; clc;

addpath '/data/tu_muellerm/utils/'

ref = trust_reference_data();

% plot colors
c_meas  = [0.8500 0.3250 0.0980];   % orange
c_model = [0      0.4470 0.7410];   % blue
c_line  = [0.9290 0.6940 0.1250];   % yellow


% =========================================================================
% 1. Proportionality of p_eff with V
% =========================================================================
r30 = ref.peff_30 ./ ref.V_30; % ref at 30V
r70 = ref.peff_70 ./ ref.V_70; % ref at 70V

in30 = ref.V_30 >= ref.V_valid(1) & ref.V_30 <= ref.V_valid(2);
in70 = ref.V_70 >= ref.V_valid(1) & ref.V_70 <= ref.V_valid(2);

fprintf('1. Proportionality of p_eff with V (%g-%g V)\n', ref.V_valid);
fprintf('   30 mm : %.1f kPa/V, sd %.1f %%\n', ...
        mean(r30(in30))/1e3, 100*std(r30(in30))/mean(r30(in30)));
fprintf('   70 mm : %.1f kPa/V, sd %.1f %%\n\n', ...
        mean(r70(in70))/1e3, 100*std(r70(in70))/mean(r70(in70)));


% =========================================================================
% 2. Depth independence of the normalised voltage response
% =========================================================================
s30 = ref.peff_30 / interp1(ref.V_30, ref.peff_30, ref.V_ref);
s70 = ref.peff_70 / interp1(ref.V_70, ref.peff_70, ref.V_ref);

V_common = intersect(ref.V_30, ref.V_70);
s30_c = interp1(ref.V_30, s30, V_common);
s70_c = interp1(ref.V_70, s70, V_common);

fprintf('2. Normalised voltage response\n');
fprintf('   %4s  %10s  %10s  %10s\n', 'V', 'res30 [%]', 'res70 [%]', 'diff [%]');
for i = 1:numel(V_common)
    ideal = V_common(i) / ref.V_ref;
    fprintf('   %4g  %10.1f  %10.1f  %10.1f\n', V_common(i), ...
            100*(s30_c(i) - ideal)/ideal, ...
            100*(s70_c(i) - ideal)/ideal, ...
            100*(s30_c(i) - s70_c(i))/s70_c(i));
end
fprintf('\n');


% =========================================================================
% 3. k(z)
% =========================================================================
[~, info] = trust_p0(ref.V_ref, ref.z_hyd, ref);
k_nodes = info.k_nodes;

fprintf('3. k(z) at the measured depths\n');
fprintf('   %4s  %10s\n', 'z', 'k [kPa]');
fprintf('   %4g  %10.1f\n', [ref.z_hyd(:), k_nodes(:)/1e3].');

plateau = ref.z_hyd >= 60 & ref.z_hyd <= 90;
fprintf('\n   plateau 60-90 mm : %.1f kPa, spread %.1f %%\n', ...
        mean(k_nodes(plateau))/1e3, ...
        100*(max(k_nodes(plateau)) - min(k_nodes(plateau))) / mean(k_nodes(plateau)));
fprintf('   20 mm relative to plateau : %+.0f %%\n\n', ...
        100*(k_nodes(1)/mean(k_nodes(plateau)) - 1));


% =========================================================================
% 4. Measured and simulated depth profiles, normalised at 50 mm
% =========================================================================
G_at_hyd = interp1(ref.z_sim, ref.G, ref.z_hyd, 'pchip');

i50       = find(ref.z_hyd == 50);
prof_meas = ref.peff_hyd / ref.peff_hyd(i50);
prof_sim  = G_at_hyd     / G_at_hyd(i50);

fprintf('4. Depth profiles normalised at 50 mm\n');
fprintf('   %4s  %9s  %9s  %9s\n', 'z', 'measured', 'model', 'diff [%]');
fprintf('   %4g  %9.3f  %9.3f  %9.1f\n', ...
        [ref.z_hyd(:), prof_meas(:), prof_sim(:), ...
         100*(prof_sim(:) - prof_meas(:))./prof_meas(:)].');
fprintf('\n');


% =========================================================================
% 5. Prediction of the voltage sweeps from the depth sweep alone
% =========================================================================

% p0 = s(V) k(z) and p_eff_pred = p0 G_sim(z), so G_sim cancels and the
% prediction reduces to (V/V_ref) times the depth-sweep value. This therefore
% tests the factorisation and the mutual consistency of the three tables. It
% does NOT test the simulation.

fprintf('5. Voltage sweeps predicted from the depth sweep\n');
fprintf('   %4s  %11s  %11s\n', 'V', 'err30 [%]', 'err70 [%]');

t30 = ref.peff_hyd(ref.z_hyd == 30);
t70 = ref.peff_hyd(ref.z_hyd == 70);

V_test = V_common(V_common >= ref.V_valid(1) & V_common <= ref.V_valid(2));
err30 = zeros(size(V_test));
err70 = zeros(size(V_test));

for i = 1:numel(V_test)
    V = V_test(i);
    pred30 = (V/ref.V_ref) * t30;
    pred70 = (V/ref.V_ref) * t70;
    meas30 = ref.peff_30(ref.V_30 == V);
    meas70 = ref.peff_70(ref.V_70 == V);
    err30(i) = 100*(pred30 - meas30)/meas30;
    err70(i) = 100*(pred70 - meas70)/meas70;
    fprintf('   %4g  %11.1f  %11.1f\n', V, err30(i), err70(i));
end

fprintf('\n   worst residual : %.1f %%\n\n', max(abs([err30 err70])));


% =========================================================================
% FIGURES
% =========================================================================

% Figure 1 (Methods): voltage scaling
figure; hold on; box on; grid on;

plot([0 90], [0 90]/ref.V_ref, 'k--', 'LineWidth', 1.0, ...
     'DisplayName', 's(V) = V/V_{ref}');
plot(ref.V_30, s30, 'o-', 'Color', c_meas, 'MarkerFaceColor', 'w', ...
     'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', '30 mm');
plot(ref.V_70, s70, 's-', 'Color', c_line, 'MarkerFaceColor', 'w', ...
     'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', '70 mm');

xlabel('Excitation voltage [V]');
ylabel('s(V) = p_{eff}(V) / p_{eff}(V_{ref})');
legend('Location', 'northwest');
set(gca, 'FontSize', 12);


% Figure 2: depth profiles normalised at 50 mm
figure; hold on; box on; grid on;

plot(ref.z_hyd, prof_meas, 'o-', 'Color', c_meas, 'MarkerFaceColor', 'w', ...
     'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Measured (hydrophone)');
plot(ref.z_hyd, prof_sim, 's--', 'Color', c_model, 'MarkerFaceColor', 'w', ...
     'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Simulated');

xlabel('Commanded focal depth [mm]');
ylabel('Normalised focal pressure (1 at 50 mm)');
legend('Location', 'northeast');
set(gca, 'FontSize', 12);


% Figure 3 (Results): required source amplitude vs commanded depth
V_plot = ref.V_ref;          % drive voltage the figure is drawn for [V]

figure; hold on; box on; grid on;

z_fine = linspace(ref.z_valid(1), ref.z_valid(2), 400);
p0_fine  = trust_p0(V_plot, z_fine,    ref);
p0_nodes = trust_p0(V_plot, ref.z_hyd, ref);

p0 = trust_p0(40, z_fine,    ref);

plot(z_fine, p0_fine/1e3, '-', 'Color', c_line, 'LineWidth', 1.8, ...
     'DisplayName', 'Interpolation, 30 V');
plot(z_fine, p0/1e3, '-', 'Color', 'r', 'LineWidth', 1.8, ...
     'DisplayName', 'Interpolation, 40 V');
plot(ref.z_hyd, p0_nodes/1e3, 'o', 'Color', c_meas, 'MarkerFaceColor', 'w', ...
     'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Calibration points');

xlabel('Commanded focal depth [mm]');
ylabel(sprintf('Required source amplitude p_0 [kPa] at %g V', V_plot));
legend('Location', 'northeast');
set(gca, 'FontSize', 12);

ref = trust_reference_data();

% Valid calibration domain
V_values = linspace(20,80,100);      % [V]
z_values = linspace(20,120,150);     % [mm]

[V_grid,z_grid] = meshgrid(V_values,z_values);

% Evaluate the actual calibration function
p0_grid = trust_p0(V_grid,z_grid,ref);   % [Pa]
