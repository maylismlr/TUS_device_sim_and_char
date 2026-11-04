%% Code to extract metrics before and after PAC
% extract_PAC_metrics.m
%
% From the data obtained after running tFUS_3D_PAC.m, obtain the tables
% found in the thesis.
%
% Compare:
%   PAC      vs free-field reference
%   No PAC   vs free-field reference
%
% Metrics:
%   - p_target
%   - p_peak
%   - peak position
%   - target-to-peak error
%   - z error
%   - -6 dB focal volume
%   - -6 dB Dice coefficient
%   - relative pressure error vs reference
%   - relative focal-volume error vs reference
%
% Author: Maylis Mueller

addpath '/data/tu_muellerm/utils/'
addpath '/data/hu_muellerm/MATLAB Add-Ons/Collections/k-Wave/k-Wave'

% ========================================================================
% LOAD FILES
% =========================================================================

% name = 'MR_HR_a2862';
name = 'muCT';

pac_steered_path = sprintf('/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/pac/%s_p_max_-10mm_45mm_PAC.mat', name);
pac_unsteered_path = sprintf('/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/pac/%s_p_max_0mm_45mm_PAC.mat', name);

muCT_steered_path = sprintf('/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/steered/%s_power2_p_max_-10mm_45mm.mat', name);
muCT_unsteered_path = sprintf('/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/unsteered/%s_power2_p_max_0mm_45mm.mat', name);

ref_steered_path = '/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/steered/ref_p_max_-10-45mm.mat';
ref_unsteered_path = '/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/unsteered/ref_p_max_0-45mm.mat';



% ========================================================================
% DEFINE CASES
% =========================================================================

cases(1).name      = 'PAC steered';
cases(1).data      = load(pac_steered_path);
cases(1).target_mm = [-10, 0, 45];
cases(1).ref_idx   = 5;

cases(2).name      = 'PAC unsteered';
cases(2).data      = load(pac_unsteered_path);
cases(2).target_mm = [0, 0, 45];
cases(2).ref_idx   = 6;

cases(3).name      = 'No PAC steered';
cases(3).data      = load(muCT_steered_path);
cases(3).target_mm = [-10, 0, 45];
cases(3).ref_idx   = 5;

cases(4).name      = 'No PAC unsteered';
cases(4).data      = load(muCT_unsteered_path);
cases(4).target_mm = [0, 0, 45];
cases(4).ref_idx   = 6;

cases(5).name      = 'Reference steered';
cases(5).data      = load(ref_steered_path);
cases(5).target_mm = [-10, 0, 45];
cases(5).ref_idx   = 5;

cases(6).name      = 'Reference unsteered';
cases(6).data      = load(ref_unsteered_path);
cases(6).target_mm = [0, 0, 45];
cases(6).ref_idx   = 6;

% ========================================================================
% OUTPUT SETTINGS
% =========================================================================

output_folder = sprintf( ...
    '/data/tu_muellerm_cloud/skulls_FUN/results/PAC_metrics/%s', name);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Figures are written to a subfolder so the CSV/MAT outputs stay easy to find.
% The .fig is kept alongside the PNG so axes limits, colour scales and fonts
% can be adjusted later for the thesis without re-running the analysis.
fig_folder = fullfile(output_folder, 'figures');

if ~exist(fig_folder, 'dir')
    mkdir(fig_folder);
end

% Save EVERYTHING printed to the command window
diary_file = fullfile(output_folder, ...
    sprintf('PAC_metrics_%s_console.txt', name));

diary(diary_file);
diary on;

fprintf('PAC analysis: %s\n', name);
fprintf('Date: %s\n\n', datetime("now"));


% ========================================================================
% GRID
% =========================================================================

dx = 0.25e-3;
dy = dx;
dz = dx;

grid_size_x = 56e-3;
grid_size_y = 56e-3;
grid_size_z = 70e-3;

Nx = roundEven(grid_size_x / dx);
Ny = roundEven(grid_size_y / dy);
Nz = roundEven(grid_size_z / dz);

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% User-coordinate depth axis:
% z = 0 mm corresponds to transducer/grid front face
z_user_mm = 1e3 * (kgrid.z_vec - kgrid.z_vec(1));

x_mm = 1e3 * kgrid.x_vec;
y_mm = 1e3 * kgrid.y_vec;


% ========================================================================
% ANALYSIS SETTINGS
% =========================================================================

% Search for the focal peak only near the intended target.
%
% This prevents a strong pressure value at the skull/interface from being
% mistaken for the therapeutic focal peak.
roi_halfwidth_mm = [15, 15, 20];   % [x y z]

num_cases = length(cases);

% Storage
p_target_MPa        = nan(num_cases,1);
p_peak_MPa          = nan(num_cases,1);

x_peak_mm           = nan(num_cases,1);
y_peak_mm           = nan(num_cases,1);
z_peak_mm           = nan(num_cases,1);

target_error_mm      = nan(num_cases,1);
z_target_error_mm    = nan(num_cases,1);

V6_mm3               = nan(num_cases,1);

mask6 = cell(num_cases,1);
peak_indices = nan(num_cases,3);


% ========================================================================
% EXTRACT BASIC METRICS FOR EACH FIELD
% =========================================================================

for i = 1:num_cases

    case_name = cases(i).name;
    p_field   = cases(i).data.p_max;

    target_mm = cases(i).target_mm;

    % ---------------------------------------------------------------------
    % Sanity check
    % ---------------------------------------------------------------------

    if ~isequal(size(p_field), [Nx Ny Nz])
        error(['Size mismatch for "%s".\n' ...
               'Pressure field = [%d %d %d]\n' ...
               'Expected       = [%d %d %d]'], ...
               case_name, ...
               size(p_field,1), size(p_field,2), size(p_field,3), ...
               Nx, Ny, Nz);
    end

    % --------------------------------------------------------------------- 
    % TARGET POSITION 
    % --------------------------------------------------------------------- 
    target_point_user = target_mm * 1e-3; 
    
    
    % Convert target z to grid coordinates 
    target_point_grid = [target_point_user(1), target_point_user(2), ...
        kgrid.z_vec(1) + target_point_user(3)];
    

    % ---------------------------------------------------------------------
    % TARGET INDICES
    % ---------------------------------------------------------------------

    [~, ix_target] = min(abs(x_mm      - target_mm(1)));
    [~, iy_target] = min(abs(y_mm      - target_mm(2)));
    [~, iz_target] = min(abs(z_user_mm - target_mm(3)));


    % ---------------------------------------------------------------------
    % PRESSURE AT INTENDED TARGET
    % ---------------------------------------------------------------------

    p_target = p_field(ix_target, iy_target, iz_target);
    p_target_MPa(i) = p_target * 1e-6;

    % --------------------------------------------------------------------- 
    % PLOT X-Z SLICE THROUGH TARGET Y 
    % --------------------------------------------------------------------- 
    prop_vec_shifted = 1e3 * (kgrid.z_vec - kgrid.z_vec(1));

    fig_field = figure;
    
    imagesc( ...
        1e3 * kgrid.x_vec, ...
        prop_vec_shifted, ...
        (1e-6 * squeeze(p_field(:, iy_target, :))).');
    
    xlabel('x-position [mm]');
    ylabel('z-position [mm]');
    axis image;
    
    title(sprintf( ...
        '%s | Target = [%.0f %.0f %.0f] mm', ...
        case_name, ...
        target_point_user(1)*1e3, ...
        target_point_user(2)*1e3, ...
        target_point_user(3)*1e3));
    
    cb = colorbar;
    title(cb, '[MPa]');

    % Save the field slice.
    fig_base = fullfile(fig_folder, ...
        sprintf('%s_%s_pressured_field', name, strrep(case_name, ' ', '_')));

    savefig(fig_field, [fig_base '.fig']);
    exportgraphics(fig_field, [fig_base '.png'], 'Resolution', 300);

    % ---------------------------------------------------------------------
    % SEARCH FOR FOCAL PEAK INSIDE TARGET-CENTRED ROI
    % ---------------------------------------------------------------------

    ix_roi = find(abs(x_mm      - target_mm(1)) <= roi_halfwidth_mm(1));
    iy_roi = find(abs(y_mm      - target_mm(2)) <= roi_halfwidth_mm(2));
    iz_roi = find(abs(z_user_mm - target_mm(3)) <= roi_halfwidth_mm(3));

    p_roi = p_field(ix_roi, iy_roi, iz_roi);

    [p_peak, linear_idx] = max(p_roi(:));

    [ix_local, iy_local, iz_local] = ...
        ind2sub(size(p_roi), linear_idx);

    ix_peak = ix_roi(ix_local);
    iy_peak = iy_roi(iy_local);
    iz_peak = iz_roi(iz_local);

    peak_indices(i,:) = [ix_peak iy_peak iz_peak];


    % Physical peak position
    x_peak_mm(i) = x_mm(ix_peak);
    y_peak_mm(i) = y_mm(iy_peak);
    z_peak_mm(i) = z_user_mm(iz_peak);

    p_peak_MPa(i) = p_peak * 1e-6;


    % ---------------------------------------------------------------------
    % FOCAL POSITION ERROR RELATIVE TO REQUESTED TARGET
    % ---------------------------------------------------------------------

    peak_position = [ ...
        x_peak_mm(i), ...
        y_peak_mm(i), ...
        z_peak_mm(i)];

    target_error_mm(i) = norm(peak_position - target_mm);

    z_target_error_mm(i) = ...
        z_peak_mm(i) - target_mm(3);


    % ---------------------------------------------------------------------
    % -6 dB FOCAL VOLUME
    % ---------------------------------------------------------------------

    vol = focal_vol_dB( ...
        p_field, ...
        dx, dy, dz, ...
        true, ...
        [ix_peak iy_peak iz_peak], ...
        'LevelDB', -6, ...
        'DBConvention', 'half');

    % Use marching-cubes volume
    V6_mm3(i) = vol.V_isosurf_mm3;


    % ---------------------------------------------------------------------
    % -6 dB CONNECTED FOCAL MASK
    % ---------------------------------------------------------------------

    mask6{i} = get_focal_mask( ...
        p_field, ...
        [ix_peak iy_peak iz_peak], ...
        0.5);


    % ---------------------------------------------------------------------
    % DISPLAY
    % ---------------------------------------------------------------------

    fprintf('\n%s\n', case_name);
    fprintf('---------------------------------------\n');
    fprintf('Target                : [%6.2f %6.2f %6.2f] mm\n', ...
        target_mm);
    fprintf('Peak                  : [%6.2f %6.2f %6.2f] mm\n', ...
        peak_position);
    fprintf('Target pressure       : %.4f MPa\n', ...
        p_target_MPa(i));
    fprintf('Peak pressure         : %.4f MPa\n', ...
        p_peak_MPa(i));
    fprintf('Target-to-peak error  : %.3f mm\n', ...
        target_error_mm(i));
    fprintf('z target error        : %+.3f mm\n', ...
        z_target_error_mm(i));
    fprintf('-6 dB focal volume    : %.3f mm^3\n', ...
        V6_mm3(i));

end


% ========================================================================
% DIFFERENCE BETWEEN PAC AND NO PAC
% =========================================================================

pair_names = {'Steered', 'Unsteered'};

pac_idx   = [1 2];
nopac_idx = [3 4];

for j = 1:2

    p_pac   = cases(pac_idx(j)).data.p_max;
    p_nopac = cases(nopac_idx(j)).data.p_max;

    target_mm = cases(pac_idx(j)).target_mm;

    [~, iy_target] = min(abs(y_mm - target_mm(2)));

    % PAC - No PAC
    diff_field = (p_pac - p_nopac) * 1e-6;   % MPa

    % X-Z slice through target y
    diff_slice = squeeze(diff_field(:, iy_target, :)).';

    fig_diff = figure;

    imagesc(x_mm, z_user_mm, diff_slice);

    xlabel('x-position [mm]', 'FontSize', 16);
    ylabel('z-position [mm]', 'FontSize', 16);

    axis image;

    title(sprintf('PAC - No PAC | %s', pair_names{j}), ...
        'FontSize', 17);

    cb = colorbar;
    title(cb, '\Delta p [MPa]');
    cb.FontSize = 14;

    set(gca, 'FontSize', 14);

    % Symmetric colour limits around zero
    cmax = max(abs(diff_slice(:)));
    clim([-cmax cmax]);

    % Blue -> white -> red
    n = 256;

    cmap = [ ...
        linspace(0,1,n/2)' linspace(0,1,n/2)' ones(n/2,1);
        ones(n/2,1) linspace(1,0,n/2)' linspace(1,0,n/2)' ];

    colormap(fig_diff, cmap);
        
    fig_base = fullfile(fig_folder, ...
        sprintf('PAC_diff_%s_%s', name, pair_names{j}));
    savefig(fig_diff, [fig_base '.fig']);
    exportgraphics(fig_diff, [fig_base '.png'], 'Resolution', 300);

end


% ========================================================================
% COMPARE EACH CASE WITH ITS FREE-FIELD REFERENCE
% =========================================================================

dice6              = nan(num_cases,1);

epsilon_p_peak_pct = nan(num_cases,1);
epsilon_p_targ_pct = nan(num_cases,1);

epsilon_V6_pct     = nan(num_cases,1);

delta_z_ref_mm     = nan(num_cases,1);
delta_r_ref_mm     = nan(num_cases,1);


for i = 1:num_cases

    r = cases(i).ref_idx;


    % ---------------------------------------------------------------------
    % DICE OF -6 dB FOCAL REGION
    % ---------------------------------------------------------------------

    A = mask6{i};
    B = mask6{r};

    dice6(i) = get_DSC(A, B);


    % ---------------------------------------------------------------------
    % PRESSURE ERROR RELATIVE TO FREE FIELD
    %
    % Signed:
    %   negative = lower than free field
    %   positive = higher than free field
    % ---------------------------------------------------------------------

    epsilon_p_peak_pct(i) = ...
        100 * ...
        (p_peak_MPa(i) - p_peak_MPa(r)) / ...
        p_peak_MPa(r);

    epsilon_p_targ_pct(i) = ...
        100 * ...
        (p_target_MPa(i) - p_target_MPa(r)) / ...
        p_target_MPa(r);


    % ---------------------------------------------------------------------
    % -6 dB FOCAL VOLUME ERROR
    % ---------------------------------------------------------------------

    epsilon_V6_pct(i) = ...
        100 * ...
        (V6_mm3(i) - V6_mm3(r)) / ...
        V6_mm3(r);


    % ---------------------------------------------------------------------
    % FOCAL POSITION RELATIVE TO FREE-FIELD PEAK
    % ---------------------------------------------------------------------

    delta_z_ref_mm(i) = ...
        z_peak_mm(i) - z_peak_mm(r);

    delta_r_ref_mm(i) = norm( ...
        [x_peak_mm(i), y_peak_mm(i), z_peak_mm(i)] - ...
        [x_peak_mm(r), y_peak_mm(r), z_peak_mm(r)]);

end


% ========================================================================
% RESULTS TABLE
% =========================================================================

Case = string({cases.name}).';

results = table( ...
    Case, ...
    p_target_MPa, ...
    p_peak_MPa, ...
    x_peak_mm, ...
    y_peak_mm, ...
    z_peak_mm, ...
    target_error_mm, ...
    z_target_error_mm, ...
    V6_mm3, ...
    dice6, ...
    epsilon_p_targ_pct, ...
    epsilon_p_peak_pct, ...
    epsilon_V6_pct, ...
    delta_z_ref_mm, ...
    delta_r_ref_mm);

disp(results);

% Save full results table
results_csv = fullfile(output_folder, ...
    sprintf('PAC_metrics_%s_results.csv', name));

writetable(results, results_csv);


% ========================================================================
% PAC IMPROVEMENT: STEERED
% =========================================================================
% The previous metric, (p_PAC - p_noPAC) / (p_ref - p_noPAC), normalised the
% PAC gain by the TOTAL target deficit relative to free field. That deficit
% contains two physically distinct contributions: phase aberration, which a
% phase-only array can correct, and transmission loss (reflection, scattering,
% absorption), which it cannot. Consequences:
%   - a perfect phase-only correction cannot reach 1, so the metric is not a
%     "fraction of the achievable" despite reading like one;
%   - it is unbounded, and changes sign whenever the skull field happens to
%     exceed free field at the target voxel;
%   - it compresses realistic corrections into the first few percent, because
%     the denominator is dominated by attenuation.
% It is therefore dropped in favour of four bounded ratios.
%
%   T_noPAC : p_t,noPAC / p_t,ref    target transmission with geometric delays
%   T_PAC   : p_t,PAC   / p_t,ref    residual target deficit after correction
%   G_t     : p_t,PAC   / p_t,noPAC  target-pressure gain      (lit.: Ptc/Pt0)
%   G_p     : p_p,PAC   / p_p,noPAC  peak-pressure gain        (lit.: Ppc/Pp0)
%
% The T_* pair makes the uncorrectable attenuation floor explicit instead of
% hiding it inside a denominator: T_PAC states how far from free field the
% target still is once the phases are right. G_t is the figure the correction
% should be judged on here, because the array is asked to place pressure at
% the commanded target, not to maximise a peak somewhere nearby. G_p is
% reported alongside only because it is the more common published benchmark;
% the two separate when the uncorrected focus is displaced, and that gap is
% itself informative.

i_pac = 1;
i_no  = 3;
i_ref = 5;

% Transmission relative to free field (bounded by 0 from below; values above
% 1 are possible and physical if the skull layer focuses at that voxel).
T_noPAC_steered = p_target_MPa(i_no)  / p_target_MPa(i_ref);
T_PAC_steered   = p_target_MPa(i_pac) / p_target_MPa(i_ref);

% Gain delivered by the correction. Denominators are the uncorrected fields,
% which are strictly positive, so these are well conditioned.
G_target_steered = p_target_MPa(i_pac) / p_target_MPa(i_no);
G_peak_steered   = p_peak_MPa(i_pac)   / p_peak_MPa(i_no);

SteeredMetric = [ ...
    "Reference target pressure [MPa]";
    "No PAC target pressure [MPa]";
    "PAC target pressure [MPa]";
    "Target transmission, no PAC [%]";
    "Target transmission, PAC [%]";
    "Target pressure gain, PAC/no PAC (Ptc/Pt0) [%]";
    "Peak pressure gain, PAC/no PAC (Ppc/Pp0) [%]";
    "No PAC Dice";
    "PAC Dice";
    "No PAC focal-position error [mm]";
    "PAC focal-position error [mm]"];

SteeredValue = [ ...
    p_target_MPa(i_ref);
    p_target_MPa(i_no);
    p_target_MPa(i_pac);
    100 * T_noPAC_steered;
    100 * T_PAC_steered;
    100 * G_target_steered;
    100 * G_peak_steered;
    dice6(i_no);
    dice6(i_pac);
    target_error_mm(i_no);
    target_error_mm(i_pac)];

steered_summary = table( ...
    SteeredMetric, ...
    SteeredValue, ...
    'VariableNames', {'Metric', 'Value'});

writetable( ...
    steered_summary, ...
    fullfile(output_folder, ...
    sprintf('PAC_metrics_%s_steered_summary.csv', name)));

% ========================================================================
% PAC IMPROVEMENT: UNSTEERED
% =========================================================================
% Same four ratios as above; see the steered block for their definition and
% for why the previous free-field-normalised recovery fraction was dropped.

i_pac = 2;
i_no  = 4;
i_ref = 6;

T_noPAC_unsteered = p_target_MPa(i_no)  / p_target_MPa(i_ref);
T_PAC_unsteered   = p_target_MPa(i_pac) / p_target_MPa(i_ref);

G_target_unsteered = p_target_MPa(i_pac) / p_target_MPa(i_no);
G_peak_unsteered   = p_peak_MPa(i_pac)   / p_peak_MPa(i_no);

UnsteeredMetric = [ ...
    "Reference target pressure [MPa]";
    "No PAC target pressure [MPa]";
    "PAC target pressure [MPa]";
    "Target transmission, no PAC [%]";
    "Target transmission, PAC [%]";
    "Target pressure gain, PAC/no PAC (Ptc/Pt0) [%]";
    "Peak pressure gain, PAC/no PAC (Ppc/Pp0) [%]";
    "No PAC Dice";
    "PAC Dice";
    "No PAC focal-position error [mm]";
    "PAC focal-position error [mm]"];

UnsteeredValue = [ ...
    p_target_MPa(i_ref);
    p_target_MPa(i_no);
    p_target_MPa(i_pac);
    100 * T_noPAC_unsteered;
    100 * T_PAC_unsteered;
    100 * G_target_unsteered;
    100 * G_peak_unsteered;
    dice6(i_no);
    dice6(i_pac);
    target_error_mm(i_no);
    target_error_mm(i_pac)];

unsteered_summary = table( ...
    UnsteeredMetric, ...
    UnsteeredValue, ...
    'VariableNames', {'Metric', 'Value'});

writetable( ...
    unsteered_summary, ...
    fullfile(output_folder, ...
    sprintf('PAC_metrics_%s_unsteered_summary.csv', name)));



% and save everything
save( ...
    fullfile(output_folder, ...
    sprintf('PAC_metrics_%s_all.mat', name)), ...
    'results', ...
    'steered_summary', ...
    'unsteered_summary', ...
    'p_target_MPa', ...
    'p_peak_MPa', ...
    'x_peak_mm', ...
    'y_peak_mm', ...
    'z_peak_mm', ...
    'target_error_mm', ...
    'z_target_error_mm', ...
    'V6_mm3', ...
    'dice6', ...
    'epsilon_p_targ_pct', ...
    'epsilon_p_peak_pct', ...
    'epsilon_V6_pct', ...
    'delta_z_ref_mm', ...
    'delta_r_ref_mm', ...
    'delta_r_ref_mm', ...
    'T_noPAC_steered', ...
    'T_PAC_steered', ...
    'G_target_steered', ...
    'G_peak_steered', ...
    'T_noPAC_unsteered', ...
    'T_PAC_unsteered', ...
    'G_target_unsteered', ...
    'G_peak_unsteered', ...
    'peak_indices');

diary off;