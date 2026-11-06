%% Plot pressure fields for all .mat files with target marker
% Script to load the pmax pressure fields computed when running
% tFUS_3D_xx.m and plotting them all with the same scale and with a cross
% at input target. Used to get Figures 3.9 and 3.10 from Master's thesis.
%
% author: Maylis Mueller

clearvars;
close all;

% =========================================================================
% SETTINGS
% =========================================================================

folder = '/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/unsteered/';
files = dir(fullfile(folder, '*.mat'));

dx = 0.25;   % [mm]
dy = dx;
dz = dx;

target = [0, 0, 45];   % [x y z] mm


% =========================================================================
% LOOP THROUGH FILES
% =========================================================================

output_folder = '/data/tu_muellerm_cloud/skulls_FUN/sensor_data/el_size_calib/steered/common_scale';

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end


% =========================================================================
% FIND GLOBAL MAXIMUM
% =========================================================================

global_max = 0;

for i = 1:numel(files)

    filepath = fullfile(folder, files(i).name);
    fname = files(i).name;

    S = load(filepath, 'p_max');

    if ~isfield(S, 'p_max')
        warning('%s does not contain p_max. Skipping.', files(i).name);
        continue;
    end


    if ~contains(fname, 'ref')
        continue
    end

    global_max = max(global_max, max(S.p_max(:)) * 1e-6);
end

fprintf('Common pressure scale: 0 to %.3f MPa\n', global_max);


for i = 1:numel(files)

    fname = files(i).name;
    filepath = fullfile(folder, fname);

    fprintf('Plotting %d/%d: %s\n', i, numel(files), fname);

    
    if contains(fname, 'muCT')
        model_name = '\muCT res. 0.125 mm';
    
    elseif contains(fname, 'CT_0.9')
        model_name = 'CT res. 0.9 mm';

    elseif contains(fname, 'ref')
        model_name = 'Reference unsteered | Target = [0 0 45] mm';
       % continue
    
    elseif contains(fname, 'CT_0.67')
        model_name = 'CT res. 0.67 mm';
    
    elseif contains(fname, 'MR_HR') && contains(fname, 'a1272')
        model_name = 'HR (0.4 mm) MRI Method 1';

    elseif contains(fname, 'MR_HR') && contains(fname, 'a2862')
        model_name = 'HR (0.4 mm) MRI Method 2';
        
    elseif contains(fname, 'MR_LR') && contains(fname, 'a1386')
        model_name = 'LR (0.75 mm) MRI Method 1';

    elseif contains(fname, 'MR_LR') && contains(fname, 'a1853')
        model_name = 'LR (0.75 mm) MRI Method 2'; 
    
    else
        [~, model_name] = fileparts(fname);
        model_name = strrep(model_name, '_', ' ');
    end

    S = load(filepath, 'p_max');

    if ~isfield(S, 'p_max')
        warning('%s does not contain p_max. Skipping.', fname);
        continue;
    end

    p_max = S.p_max;

    [Nx, Ny, Nz] = size(p_max);

    % ---------------------------------------------------------------------
    % Coordinate vectors
    % ---------------------------------------------------------------------

    % x and y are centred around zero
    x_mm = ((0:Nx-1) - (Nx-1)/2) * dx;
    y_mm = ((0:Ny-1) - (Ny-1)/2) * dy;

    % z = 0 at the transducer surface
    z_mm = (0:Nz-1) * dz;

    % Find y = 0 slice
    [~, iy] = min(abs(y_mm - target(2)));

    % x-z pressure field
    p_xz = squeeze(p_max(:, iy, :));

    % ---------------------------------------------------------------------
    % Plot
    % ---------------------------------------------------------------------

    fig = figure('Color', 'w');
    ax = axes(fig);
    
    imagesc(ax, x_mm, z_mm, (p_xz * 1e-6).');
    axis(ax, 'image');

    hold on;

    % Target [-10, 0, 40] mm
    plot(target(1), target(3), 'kx', ...
        'MarkerSize', 18, 'MarkerEdgeColor', 'r', ...
        'LineWidth', 2.5);

    xlabel('x-position [mm]');
    ylabel('z-position [mm]');

    % Same pressure scale
    clim(ax, [0 global_max]);

    cb = colorbar;
    ylabel(cb, 'Pressure [MPa]');

    title(strrep(model_name, '_', '\_'), 'FontSize', 17);

    set(gca, 'FontSize', 14);

end

% =========================================================================
% APPLY COMMON SCALE AND DISPLAY
% =========================================================================

for i = 1:numel(files)

    fprintf('Processing %d/%d: %s\n', ...
        i, numel(files), output_names{i});

    % OPEN VISIBLE
    fig = openfig(files{i}, 'visible');

    hImg = findobj(fig, 'Type', 'image');
    ax = ancestor(hImg(1), 'axes');

    % Same pressure scale
    clim(ax, [0 global_max]);

    % Typography
    ax.Title.FontSize = 17;

    cb = findobj(fig, 'Type', 'ColorBar');
    if ~isempty(cb)
        cb(1).Title.FontSize = 14;
    end

    % Add target marker to the pressure-field axes
    hold(ax, 'on');

    plot(ax, target(1), target(3), 'x', ...
        'MarkerSize', 18, ...
        'MarkerEdgeColor', 'r', ...
        'LineWidth', 2.5);

    hold(ax, 'off');

end

disp('All pressure fields exported with the same color scale.');