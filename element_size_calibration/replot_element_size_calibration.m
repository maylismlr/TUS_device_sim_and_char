%% replot_element_size_calibration.m
%
% Load axial profiles (recorded and simulated) and tables obtained after element size calibration obtained after running fit_element_size_profile.m, then plot
% and save Figures 3.2 and 3.3, as well as Supp. Figures 6.8 and 6.9.
%
% Author: Maylis Mueller

clearvars; close all; clc;
addpath '/data/tu_muellerm/utils';

% =========================================================================
% SETTINGS
% =========================================================================
% Nominal size
element_size = 'nom_size';
selected_element_size_mm = 2.58;

% Best fit size
% element_size = '39';
% selected_element_size_mm = 3.9;

if strcmp(element_size, '39')
    mat_files = {
        '/data/tu_muellerm/element_size_calibration/element_size_calibration_figures_39/afp_element_size_fit_hydrophone_broad.mat'
        '/data/tu_muellerm/element_size_calibration/element_size_calibration_figures_39/afp_element_size_fit_hydrophone_narrow.mat'
    };
elseif strcmp(element_size, 'nom_size')
    mat_files = {
        '/data/tu_muellerm/element_size_calibration/element_size_calibration_figures_nom_size/afp_element_size_fit_hydrophone.mat',...
        };
else
    error('Element size should be either "39" or "nom_size"')
end


save_figures = false;
figure_dir = sprintf('element_size_calibration_figures_%s', element_size);

if save_figures && ~exist(figure_dir,'dir') 
    mkdir(figure_dir); 
end

% colors
orange = [0.8500 0.3250 0.0980];
blue   = [0 0.4470 0.7410];


% =========================================================================
% LOAD + MERGE SWEEPS
% =========================================================================
sizes = [];
sim_peak_p = [];
sim_peak_z = [];
sim_profiles = [];
l2 = [];
mean_l2 = [];
mean_pos = [];
mean_p_err = [];

for k = 1:numel(mat_files)

    D = load(mat_files{k});

    if k == 1
        hydro = D.hydro;
        focus_values_mm = D.focus_values_mm;
        z_vec = D.z_vec;
        element_pitch = D.element_pitch;
        element_size_nominal = D.element_size_nominal;

        if isfield(D,'pressure_reference_focus_mm')
            pressure_reference_focus_mm = D.pressure_reference_focus_mm;
        else
            pressure_reference_focus_mm = 30;
        end
    end

    sizes        = [sizes; D.element_size_values_mm(:)];
    sim_peak_p   = cat(1,sim_peak_p,D.simulated_peak_pressure);
    sim_peak_z   = cat(1,sim_peak_z,D.simulated_peak_z);
    sim_profiles = cat(1,sim_profiles,D.simulated_axial_pressure);
    l2           = cat(1,l2,D.l2_error);

    mean_l2   = [mean_l2; D.mean_l2_error(:)];
    mean_pos  = [mean_pos; D.mean_focal_position_error_mm(:)];
    mean_p_err = [mean_p_err; D.mean_peak_pressure_error_pct(:)];

end


% =========================================================================
% SORT + REMOVE DUPLICATES
% =========================================================================
[sizes,idx] = sort(sizes);

sim_peak_p   = sim_peak_p(idx,:);
sim_peak_z   = sim_peak_z(idx,:);
sim_profiles = sim_profiles(idx,:,:);
l2           = l2(idx,:);
mean_l2      = mean_l2(idx);
mean_pos     = mean_pos(idx);
mean_p_err   = mean_p_err(idx);

[~,idx] = unique(round(sizes,6),'stable');

sizes        = sizes(idx);
sim_peak_p   = sim_peak_p(idx,:);
sim_peak_z   = sim_peak_z(idx,:);
sim_profiles = sim_profiles(idx,:,:);
l2           = l2(idx,:);
mean_l2      = mean_l2(idx);
mean_pos     = mean_pos(idx);
mean_p_err   = mean_p_err(idx);


% =========================================================================
% SELECT size
% =========================================================================
[d,selected_ind] = min(abs(sizes-selected_element_size_mm));

if d > 1e-6
    error('%.3f mm not found. Available sizes: %s', ...
        selected_element_size_mm, mat2str(sizes'));
end


% =========================================================================
% HYDROPHONE -> FOCUS INDEX
% =========================================================================
n_scans = numel(hydro);
n_focus = numel(focus_values_mm);

hydro_focus_ind = zeros(n_scans,1);

for k = 1:n_scans
    [~,hydro_focus_ind(k)] = ...
        min(abs(focus_values_mm-hydro(k).focus_mm));
end


% =========================================================================
% NORMALIZED PEAK PRESSURE
% =========================================================================
[~,ref_ind] = ...
    min(abs(focus_values_mm-pressure_reference_focus_mm));

hydro_peak = nan(1,n_focus);

for k = 1:n_focus
    mask = hydro_focus_ind == k;
    hydro_peak(k) = mean([hydro(mask).p_peak],'omitnan');
end

hydro_peak_norm = hydro_peak / hydro_peak(ref_ind);
sim_peak_norm   = sim_peak_p ./ sim_peak_p(:,ref_ind);


% =========================================================================
% HYDROPHONE FOCAL POSITIONS
% =========================================================================
hydro_peak_z_mm = nan(1,n_focus);

for k = 1:n_focus
    mask = hydro_focus_ind == k;
    hydro_peak_z_mm(k) = ...
        mean([hydro(mask).z_peak],'omitnan')*1e3;
end


% =========================================================================
% 1. L2 ERROR
% =========================================================================
fig = figure;
plot(sizes,mean_l2,'o-', ...
    'Color',blue, ...
    'LineWidth',1.5, ...
    'MarkerSize',7, ...
    'MarkerFaceColor',blue);
hold on;

plot(sizes(selected_ind),mean_l2(selected_ind),'x', ...
    'Color',orange, ...
    'MarkerSize',15,'LineWidth',2);

xline(element_pitch*1e3,'k:','Pitch');

grid on;
xlabel('Effective element size (mm)', 'FontSize',14);
ylabel('Mean relative L_2 error', 'FontSize',14);
title('Hydrophone profile-shape error', 'FontSize',14);

if save_figures
    saveplot(fig,'element_size_L2_error');
end


% =========================================================================
% 2. PEAK-PRESSURE ERROR
% =========================================================================
fig = figure;
plot(sizes,mean_p_err,'o-', ...
    'Color',blue, ...
    'LineWidth',1.5, ...
    'MarkerSize',7, ...
    'MarkerFaceColor',blue);
hold on;

plot(sizes(selected_ind),mean_p_err(selected_ind),'x', ...
    'Color',orange, ...
    'MarkerSize',15,'LineWidth',2);

grid on;
xlabel('Effective element size (mm)', 'FontSize',14);
ylabel('Mean normalized peak-pressure error (\epsilon_p) (%)', 'FontSize',14);
title(sprintf('Peak-pressure depth dependence (normalized at %.0f mm)', ...
    pressure_reference_focus_mm), 'FontSize',14);

if save_figures
    saveplot(fig,'element_size_peak_pressure_error');
end


% =========================================================================
% 3. FOCAL-POSITION ERROR
% =========================================================================
fig = figure;
plot(sizes,mean_pos,'o-', ...
    'Color',blue, ...
    'LineWidth',1.5, ...
    'MarkerSize',7, ...
    'MarkerFaceColor',blue);
hold on;

plot(sizes(selected_ind),mean_pos(selected_ind),'x', ...
    'Color',orange, ...
    'MarkerSize',15,'LineWidth',2);

grid on;
xlabel('Effective element size (mm)', 'FontSize',14);
ylabel('Mean axial focal-position error (\epsilon_z) (mm)', 'FontSize',14);
title('Hydrophone focal-position error', 'FontSize',14);

if save_figures
    saveplot(fig,'element_size_focal_position_error');
end


% =========================================================================
% 4. NORMALIZED PEAK PRESSURE
% =========================================================================

% =========================================================================
% 4. NORMALIZED PEAK PRESSURE BEFORE INPUT PRESSURE SCALING
% =========================================================================
fig = figure;

plot(focus_values_mm, sim_peak_norm(selected_ind,:), 'o-', ...
    'Color', blue, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', blue);
hold on;

plot(focus_values_mm, hydro_peak_norm, '^-', ...
    'Color', orange, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', orange);

grid on;
box on;

xlabel('Commanded focal depth [mm]', 'FontSize', 14);
ylabel('Peak pressure [MPa] normalized at 30 mm', 'FontSize', 14);

legend('Simulation', 'Hydrophone', ...
    'Location', 'northeast', 'FontSize', 12);

title('Normalized peak pressure before input pressure scaling', ...
    'FontSize', 14);

if save_figures
    saveplot(fig, 'normalized_peak_pressure');
end


% =========================================================================
% 5. FOCAL POSITION
% =========================================================================
fig = figure;

plot(focus_values_mm,hydro_peak_z_mm,'o-', ...
    'Color',orange, ...
    'LineWidth',1.5,'MarkerSize',7, ...
    'MarkerFaceColor',orange);
hold on;

plot(focus_values_mm,sim_peak_z(selected_ind,:)*1e3,'^-', ...
    'Color',blue, ...
    'LineWidth',1.5,'MarkerSize',7, ...
    'MarkerFaceColor',blue);

plot(focus_values_mm,focus_values_mm,'k:','LineWidth',1.2);

grid on;
xlabel('Commanded focal depth (mm)', 'FontSize',14);
ylabel('Measured / simulated peak z position (mm)', 'FontSize',14);
legend('Hydrophone','Simulation','Commanded = realised', ...
    'Location','best', 'FontSize',14);
title('Focal position before delay calibration', 'FontSize',14);

if save_figures
    saveplot(fig,'focal_position_comparison');
end


% =========================================================================
% 6. AXIAL PROFILE COMPARISONS AT SELECTED SIZE
% =========================================================================
fwhm_h = nan(n_scans);

fwhm_s = nan(n_scans);

for scan_ind = 1:n_scans

    focus_ind = hydro_focus_ind(scan_ind);

    z_hyd = hydro(scan_ind).z(:);
    p_hyd = hydro(scan_ind).p(:);

    p_sim = squeeze(sim_profiles(selected_ind,focus_ind,:));

    p_sim = interp1(z_vec(:),p_sim,z_hyd,'linear',NaN);

    valid = isfinite(p_hyd) & isfinite(p_sim);

    z = z_hyd(valid)*1e3;      % mm
    p_h = p_hyd(valid);
    p_s = p_sim(valid);

    % Normalize profiles
    p_h = p_h/max(p_h);
    p_s = p_s/max(p_s);

    % -------------------------------------------------------------
    % FWHM
    % -------------------------------------------------------------
    % Find peak indices
    [~, peak_h] = max(p_h);
    [~, peak_s] = max(p_s);
    
    % FWHM around each peak
    [fwhm_h(scan_ind), zL_h, zR_h] = get_fwhm(z, p_h, peak_h);
    [fwhm_s(scan_ind), zL_s, zR_s] = get_fwhm(z, p_s, peak_s);

    % -------------------------------------------------------------
    % Plot
    % -------------------------------------------------------------
    fig = figure;

    % Profiles
    plot(z,p_h,'-', ...
        'Color',orange,'LineWidth',1.5);
    hold on;
    
    plot(z,p_s,'--', ...
        'Color',blue,'LineWidth',1.5);
    
    xline(hydro(scan_ind).focus_mm,'k:','Commanded');
    
    % Hydrophone FWHM
    plot([zL_h zR_h],[0.5 0.5],'-', ...
        'Color',orange,'LineWidth',3, ...
        'HandleVisibility','off');
    
    plot([zL_h zL_h],[0 0.5],':', ...
        'Color',orange,'LineWidth',1.2, ...
        'HandleVisibility','off');
    
    plot([zR_h zR_h],[0 0.5],':', ...
        'Color',orange,'LineWidth',1.2, ...
        'HandleVisibility','off');
    
    plot([zL_h zR_h],[0.5 0.5],'o', ...
        'Color',orange, ...
        'MarkerFaceColor',orange, ...
        'HandleVisibility','off');
    
    % Simulation FWHM
    plot([zL_s zR_s],[0.5 0.5],'--', ...
        'Color',blue,'LineWidth',3, ...
        'HandleVisibility','off');
    
    plot([zL_s zL_s],[0 0.5],':', ...
        'Color',blue,'LineWidth',1.2, ...
        'HandleVisibility','off');
    
    plot([zR_s zR_s],[0 0.5],':', ...
        'Color',blue,'LineWidth',1.2, ...
        'HandleVisibility','off');
    
    plot([zL_s zR_s],[0.5 0.5],'o', ...
        'Color',blue, ...
        'MarkerFaceColor',blue, ...
        'HandleVisibility','off');
    
    grid on;
    xlabel('z position (mm)', 'FontSize',14);
    ylabel('Normalized pressure', 'FontSize',14);
    ylim([0 1.05]);
    
    legend('Hydrophone','Simulation','Location','best', 'FontSize',14);
    
    title(sprintf(['Focus %.0f mm: relative L_2 = %.4f'], ...
        hydro(scan_ind).focus_mm, ...
        l2(selected_ind,scan_ind)), 'FontSize',14);
    
    if save_figures
        saveplot(fig,sprintf('profile_focus_%gmm_scan_%02d',...
            hydro(scan_ind).focus_mm,scan_ind));
    end
end


% =========================================================================
% LOCAL FUNCTION
% =========================================================================
function saveplot(fig,name)

    save_figures = evalin('base','save_figures');
    if ~save_figures, return; end

    folder = evalin('base','figure_dir');

    savefig(fig,fullfile(folder,[name '.fig']));
    exportgraphics(fig,fullfile(folder,[name '.png']),...
       'Resolution',300);

end