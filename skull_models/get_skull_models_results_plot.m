%% Batch process all .mat files in a folder and compare -6 dB focal volumes
% Script to extract and plot comparison metrics across configurations / skull 
% models. Used to create Figure 3.11 and Table 3.4 from Master's thesis.
% 
% NB: here Methode 1 is "lin", Methode 2 is "hist"
%
% author: Maylis Mueller

clearvars;
addpath '/data/tu_muellerm/utils/'

folder = '/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/steered/';
files = dir(fullfile(folder, '*.mat'));

output_folder = '/data/tu_muellerm/skull_models/skull_model_results/test/';

if isempty(files)
    error('No .mat files found in %s', folder);
end

nFiles = numel(files);

results = struct('filename', {}, ...
                 'V_focal_mm3', {}, ...
                 'V_isosurf_mm3', {}, ...
                 'N_voxels', {}, ...
                 'p_max_peak', {}, ...
                 'x_peak_mm', {}, ...
                 'y_peak_mm', {}, ...
                 'z_peak_mm', {});

% Store -6 dB focal masks
mask6 = cell(nFiles,1);

% set dx so that it corresponds to dx of grid when file was created
dx = 0.25e-3;
dy = dx;
dz = dx;


% =========================================================================
% PROCESS FILES
% =========================================================================

for i = 1:nFiles

    fpath = fullfile(files(i).folder, files(i).name);
    fprintf('\nProcessing file %d/%d: %s\n', i, nFiles, files(i).name);

    S = load(fpath);
    p_max_3D = S.p_max;

    % Focal volume
    vol = focal_vol_dB(p_max_3D, dx, dy, dz, false);

    % Find focal peak using same z-search as focal_vol_dB
    sub_vol = p_max_3D(:,:,80:end);
    [~, max_idx] = max(sub_vol(:));
    [ix_peak, iy_peak, iz_peak_sub] = ind2sub(size(sub_vol), max_idx);
    iz_peak = iz_peak_sub + 79;

    % -6 dB connected focal mask
    mask6{i} = get_focal_mask(p_max_3D, [ix_peak iy_peak iz_peak], 0.5);

    % Store results
    results(end+1).filename = files(i).name;
    results(end).V_focal_mm3 = vol.V_focal_mm3;
    results(end).V_isosurf_mm3 = vol.V_isosurf_mm3;
    results(end).N_voxels = vol.N_voxels;
    results(end).p_max_peak = max(sub_vol(:));

    x_peak_mm = (ix_peak - (size(p_max_3D,1)+1)/2) * dx * 1e3;
    y_peak_mm = (iy_peak - (size(p_max_3D,2)+1)/2) * dy * 1e3;
    z_peak_mm = (iz_peak - 1) * dz * 1e3;
    
    results(end).x_peak_mm = x_peak_mm;
    results(end).y_peak_mm = y_peak_mm;
    results(end).z_peak_mm = z_peak_mm;
    
     
end


% =========================================================================
% CONVERT RESULTS TO TABLE for easy viewing/export
% =========================================================================

resultsTable = struct2table(results);

resultsTable.Properties.VariableNames = {'FileName', ...
                                         'VoxelCount_Volume_mm3', ...
                                         'Isosurface_Volume_mm3', ...
                                         'N_Voxels', ...
                                         'p_max_peak', ...
                                         'x_peak_mm', ...
                                         'y_peak_mm', ...
                                         'z_peak_mm'};


disp(resultsTable);

% writetable(resultsTable, fullfile(output_folder, 'focal_volume_results.csv'));
% writetable(resultsTable, fullfile(output_folder, 'focal_volume_results.xlsx'));


% =========================================================================
% COMPARISON TABLE WITH PERCENTAGE DIFFERENCES RELATIVE TO muCT
% =========================================================================

volMetric = 'Isosurface_Volume_mm3';

names = resultsTable.FileName;

% Identify files
idx_MR_HR_meth1 = find(contains(names,'MR_HR') & contains(names,'1272'));
idx_MR_HR_meth2 = find(contains(names,'MR_HR') & contains(names,'2862'));

idx_MR_LR_meth1 = find(contains(names,'MR_LR') & contains(names,'1386'));
idx_MR_LR_meth2 = find(contains(names,'MR_LR') & contains(names,'1853'));

% idx_CT09 = find(contains(names,'CT_0.9'));
idx_CT067 = find(contains(names,'CT_0.67'));

idx_muCT = find(contains(names,'muCT'));
idx_free_field = find(contains(names,'ref'));


% =========================================================================
% EXTRACT VOLUMES
% =========================================================================

V_MR_HR_meth1 = resultsTable.(volMetric)(idx_MR_HR_meth1);
V_MR_HR_meth2 = resultsTable.(volMetric)(idx_MR_HR_meth2);

V_MR_LR_meth1 = resultsTable.(volMetric)(idx_MR_LR_meth1);
V_MR_LR_meth2 = resultsTable.(volMetric)(idx_MR_LR_meth2);

% V_CT09 = resultsTable.(volMetric)(idx_CT09);
V_CT067 = resultsTable.(volMetric)(idx_CT067);

V_muCT = resultsTable.(volMetric)(idx_muCT);
V_free_field = resultsTable.(volMetric)(idx_free_field);


% =========================================================================
% PERCENTAGE VOLUME DIFFERENCES RELATIVE TO muCT
% =========================================================================

pct_MR_HR_meth1 = 100*(V_MR_HR_meth1 - V_muCT)/V_muCT;
pct_MR_HR_meth2 = 100*(V_MR_HR_meth2 - V_muCT)/V_muCT;

pct_MR_LR_meth1 = 100*(V_MR_LR_meth1 - V_muCT)/V_muCT;
pct_MR_LR_meth2 = 100*(V_MR_LR_meth2 - V_muCT)/V_muCT;

% pct_CT09 = 100*(V_CT09 - V_muCT)/V_muCT;
pct_CT067 = 100*(V_CT067 - V_muCT)/V_muCT;

pct_free_field = 100*(V_free_field - V_muCT)/V_muCT;
pct_muCT = 0;


% =========================================================================
% PERCENTAGE PMAX DIFFERENCE RELATIVE TO muCT
% =========================================================================

pmax_muCT = resultsTable.p_max_peak(idx_muCT);

pct_pmax = 100 * ([
    resultsTable.p_max_peak(idx_free_field);
    resultsTable.p_max_peak(idx_muCT);
    resultsTable.p_max_peak(idx_MR_HR_meth1);
    resultsTable.p_max_peak(idx_MR_HR_meth2);
    resultsTable.p_max_peak(idx_MR_LR_meth1);
    resultsTable.p_max_peak(idx_MR_LR_meth2);
    resultsTable.p_max_peak(idx_CT067)
] - pmax_muCT) / pmax_muCT;


% =========================================================================
% TRANSMISSION RELATIVE TO FREE-FIELD
% =========================================================================

transmission_ref = resultsTable.p_max_peak(idx_free_field);

transmission = 100 * ([
    resultsTable.p_max_peak(idx_free_field);
    resultsTable.p_max_peak(idx_muCT);
    resultsTable.p_max_peak(idx_MR_HR_meth1);
    resultsTable.p_max_peak(idx_MR_HR_meth2);
    resultsTable.p_max_peak(idx_MR_LR_meth1);
    resultsTable.p_max_peak(idx_MR_LR_meth2);
    resultsTable.p_max_peak(idx_CT067)
])  / transmission_ref;


% =========================================================================
% DICE SIMILARITY COEFFICIENT RELATIVE TO muCT
% =========================================================================

DSC_vs_muCT = nan(nFiles,1);

mask_muCT = mask6{idx_muCT};

for i = 1:nFiles

    if ~isequal(size(mask6{i}), size(mask_muCT))
        warning('Grid size mismatch for %s. DSC not calculated.', ...
            files(i).name);
        continue
    end

    DSC_vs_muCT(i) = get_DSC(mask6{i}, mask_muCT);
end

% Add DSC to full results table
resultsTable.DSC_vs_muCT = DSC_vs_muCT;


% =========================================================================
% AXIAL FOCAL DISPLACEMENT RELATIVE TO FREE FIELDD
% =========================================================================

x_ff = resultsTable.x_peak_mm(idx_muCT);
y_ff = resultsTable.y_peak_mm(idx_muCT);
z_ff = resultsTable.z_peak_mm(idx_muCT);

resultsTable.Epsilon_r_vs_muCT_mm = sqrt( ...
    (resultsTable.x_peak_mm - x_ff).^2 + ...
    (resultsTable.y_peak_mm - y_ff).^2 + ...
    (resultsTable.z_peak_mm - z_ff).^2);

disp(resultsTable);


% =========================================================================
% EXTRACT MR a-VALUES FOR LABELS
% =========================================================================

mr_names = names(idx_MR_HR_meth1);

a_values = regexp(mr_names,'a([\d\.]+)','tokens');
a_values = cellfun(@(c) str2double(c{1}{1}), a_values);


% =========================================================================
% BUILD COMPARISON TABLE
% =========================================================================

comparisonLabels = {
    {'Free field'}
    {'muCT 0.125mm'}
    {'MR HR 0.4mm', 'Method 1'}
    {'MR HR 0.4mm', 'Method 2'}
    {'MR LR 0.75mm', 'Method 1'}
    {'MR LR 0.75mm', 'Method 2'}
    {'CT 0.67mm'}
};

% {'CT 0.98'}

volumes = [
    V_free_field
    V_muCT
    V_MR_HR_meth1
    V_MR_HR_meth2
    V_MR_LR_meth1
    V_MR_LR_meth2
    V_CT067
];

% V_CT09

pctDiffs = [
    pct_free_field
    pct_muCT
    pct_MR_HR_meth1
    pct_MR_HR_meth2
    pct_MR_LR_meth1
    pct_MR_LR_meth2
    pct_CT067
];

% pct_CT09

diceComparison = [
    DSC_vs_muCT(idx_free_field)
    DSC_vs_muCT(idx_muCT)
    DSC_vs_muCT(idx_MR_HR_meth1)
    DSC_vs_muCT(idx_MR_HR_meth2)
    DSC_vs_muCT(idx_MR_LR_meth1)
    DSC_vs_muCT(idx_MR_LR_meth2)
    DSC_vs_muCT(idx_CT067)
];

DSC_plot = 100 * diceComparison;


epsilon_r = [
    resultsTable.Epsilon_r_vs_muCT_mm(idx_free_field)
    resultsTable.Epsilon_r_vs_muCT_mm(idx_muCT)
    resultsTable.Epsilon_r_vs_muCT_mm(idx_MR_HR_meth1)
    resultsTable.Epsilon_r_vs_muCT_mm(idx_MR_HR_meth2)
    resultsTable.Epsilon_r_vs_muCT_mm(idx_MR_LR_meth1)
    resultsTable.Epsilon_r_vs_muCT_mm(idx_MR_LR_meth2)
    resultsTable.Epsilon_r_vs_muCT_mm(idx_CT067)
];


pctDiffTable = table( ...
    cellfun(@(c) strjoin(c,' '), comparisonLabels, 'UniformOutput', false), ...
    volumes, pctDiffs, pct_pmax, diceComparison, epsilon_r, ...
    'VariableNames', {'Comparison', ...
                      'Volume_mm3', ...
                      'PctDiff_vs_muCT', ...
                      'PctPmaxDiff_vs_muCT', ...
                      'DSC_vs_muCT', ...
                      'Epsilon_r_vs_muCT_mm'});

disp(pctDiffTable);


% Optional export
% writetable(pctDiffTable, fullfile(output_folder,'focal_volume_pct_diff_muCT.csv'));


% =========================================================================
% PLOT VOLUME BARS + PERCENTAGE DIFFERENCE DOTS
% =========================================================================

fig = figure('Position', [100 100 1100 700]);

% Typography
fs_ticks  = 14;
fs_labels = 16;
fs_legend = 14;
fs_title  = 17;
fs_text   = 13;

% LEFT AXIS: focal volume

yyaxis left

bar(volumes)
ylabel('-6 dB Focal Volume (mm^3)', 'FontSize', fs_labels)

ax = gca;
ax.FontSize = fs_ticks;
ax.XTick = 1:length(volumes);
ax.XTickLabel = [];

yl = ylim;
offset = 0.08 * range(yl);

for k = 1:numel(comparisonLabels)
    text(k, yl(1)-offset, comparisonLabels{k}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'Rotation', 45, ...
        'FontSize', fs_ticks);
end

ax.Position(2) = ax.Position(2) + 0.1;
ax.Position(4) = ax.Position(4) - 0.1;

grid on
hold on


% RIGHT AXIS: differences and DSC and transmission
epsilon = char(hex2dec('03B5'));
mu =  char(hex2dec('00B5'));
yyaxis right

% focal volume difference
scatter(1:length(pctDiffs), pctDiffs, 90, 'filled');
hold on

% pmax difference
% scatter(1:length(pct_pmax), pct_pmax, 120, '*');

% DSC
scatter(1:length(DSC_plot), DSC_plot, 90, 'square', 'filled');

% Transmission
yellow  = [0.9290 0.6940 0.1250];   % yellow
scatter(1:length(transmission), transmission, 90, '*', 'MarkerEdgeColor',yellow, 'MarkerFaceColor', yellow);

ylabel(sprintf('DSC / %s_V rel. to %sCT and \nTransmission rel. to free-field (%%)', epsilon, mu), 'FontSize', fs_labels)
ylim([-40 150])
ax = gca;
ax.FontSize = fs_ticks;

yline(0, '--', 'LineWidth', 1.2);


% Value labels
for i = 1:length(pctDiffs)

    text(i, pctDiffs(i)-1, sprintf('%s_V=%.0f%%', epsilon, pctDiffs(i)), ...
        'FontWeight', 'bold', ...
        'FontSize', fs_text, ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'right');

    text(i, DSC_plot(i)+1, sprintf('DSC=%.0f%%', DSC_plot(i)), ...
        'FontWeight', 'bold', ...
        'FontSize', fs_text, ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'right');

    text(i, transmission(i), sprintf('T=%.0f%%', transmission(i)), ...
        'FontWeight', 'bold', ...
        'FontSize', fs_text, ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'right');
end


% Legend and title
lgd = legend({'Focal Volume', ... % 
              'Volume difference vs \muCT', ...
              'DSC vs \muCT', ...
              'Transmission vs free-field'}, ...
              'Location', 'northwest');

lgd.FontSize = fs_legend;

title('Comparison of Skull Models', ...
    'FontSize', fs_title, ...
    'FontWeight', 'bold');


% Export
bar_fig = gcf;

% exportgraphics(bar_fig, fullfile(output_folder, 'volume_comp.png'), ...
%     'Resolution', 300);