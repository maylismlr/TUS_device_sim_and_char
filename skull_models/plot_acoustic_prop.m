%% Plot acoustic properties of evry skull model
% Script to plot acoustic properties of evry skull model obtained after
% running the acoustical_properties_from_imaging.ipynb file. Used to obtain
% Figure 3.8 from Master's thesis.
%
% author: Maylis Mueller

clearvars;
addpath '/data/tu_muellerm/utils/'
addpath '/data/hu_muellerm/MATLAB Add-Ons/Collections/k-Wave/k-Wave'


% =========================================================================
% FOLDER CONTAINING SKULL ACOUSTIC PROPERTY FILES
% =========================================================================

skull_folder = "/data/tu_muellerm_cloud/skulls_FUN/skulls_params/final";

output_folder = fullfile(skull_folder, 'acoustic_property_plots');

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

files = dir(fullfile(skull_folder, '*.mat'));

fprintf('Found %d .mat files\n', numel(files));

dx_mm = 0.25;
dz_mm = 0.25;


% =========================================================================
% FIND GLOBAL MAXIMUM VALUES ACROSS ALL SKULL MODELS
% =========================================================================

max_density = -Inf;
max_speed   = -Inf;
max_alpha   = -Inf;

for k = 1:numel(files)

    file_path = fullfile(files(k).folder, files(k).name);

    vars = whos('-file', file_path);
    var_names = {vars.name};

    required = {'density', 'speed_of_sound', 'alpha_coeff'};

    if contains(files(k).name, 'muCT')

        S_tmp = load(file_path, 'seg');
        seg = S_tmp.seg;

        density = 1000 * ones(size(seg));
        density(seg == 1) = 2000;

        speed_of_sound = 1480 * ones(size(seg));
        speed_of_sound(seg == 1) = 3100;

        alpha_coeff = 2.17e-3 * ones(size(seg));
        alpha_coeff(seg == 1) = 2.7;

    else

        if ~all(ismember(required, var_names))
            continue;
        end

        S_tmp = load(file_path, required{:});

        density        = S_tmp.density;
        speed_of_sound = S_tmp.speed_of_sound;
        alpha_coeff    = S_tmp.alpha_coeff;

    end

    max_density = max(max_density, max(density(:)));
    max_speed   = max(max_speed, max(speed_of_sound(:)));
    max_alpha   = max(max_alpha, max(alpha_coeff(:)));

    clear S_tmp density speed_of_sound alpha_coeff
end

fprintf('Global maximum density: %.2f kg/m^3\n', max_density);
fprintf('Global maximum speed of sound: %.2f m/s\n', max_speed);
fprintf('Global maximum alpha: %.4f dB/(MHz^2 cm)\n', max_alpha);


% =========================================================================
% LOOP OVER SKULL MODELS
% =========================================================================


for k = 1:numel(files)

    file_path = fullfile(files(k).folder, files(k).name);

    fprintf('Processing %s\n', files(k).name);

    % Check what is inside before loading everything
    vars = whos('-file', file_path);
    var_names = {vars.name};

    required = {'density', 'speed_of_sound', 'alpha_coeff'};
    % 
    % if ~all(ismember(required, var_names))
    %     fprintf('  Skipping: required acoustic properties not found.\n');
    %     continue;
    % end

    if contains(files(k).name, 'mCT')
    
        fprintf('  Constant-property mCT model\n');
    
        S_tmp = load(file_path, 'seg');
        seg = S_tmp.seg;
    
        % Water background
        S.density        = 1000 * ones(size(seg));
        S.speed_of_sound = 1480 * ones(size(seg));
        S.alpha_coeff    = 2.17e-3 * ones(size(seg));
    
        % Constant bone properties
        S.density(seg == 1)        = 2000;
        S.speed_of_sound(seg == 1) = 3100;
        S.alpha_coeff(seg == 1)    = 2.7;
    
    else
    
        % Check that heterogeneous acoustic properties exist
        if ~all(ismember(required, var_names))
            fprintf('  Skipping: required acoustic properties not found.\n');
            continue;
        end
    
        S = load(file_path, required{:});
    
    end

    % Central y slice -> x-z plane
    iy = round(size(S.density, 2) / 2);
    
    density_slice = squeeze(S.density(:, iy, :));
    c_slice       = squeeze(S.speed_of_sound(:, iy, :));
    alpha_slice   = squeeze(S.alpha_coeff(:, iy, :));
    
    Nx = size(density_slice, 1);
    Nz = size(density_slice, 2);
    
    x_mm = ((0:Nx-1) - (Nx-1)/2) * dx_mm;
    z_mm = (0:Nz-1) * dz_mm;
    
    % Keep only first 25 mm in z
    z_keep = z_mm <= 25;
    
    z_mm          = z_mm(z_keep);
    density_slice = density_slice(:, z_keep);
    c_slice       = c_slice(:, z_keep);
    alpha_slice   = alpha_slice(:, z_keep);



    % Clean display name
    fname = files(k).name;
    
    if contains(fname, 'muCT')
        model_name = 'micro-CT res. 0.125 mm';
    
    elseif contains(fname, 'CT_0.9')
        model_name = 'CT res. 0.9 mm';
    
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

    % =========================================================================
    % PLOT
    % =========================================================================
    
    fig = figure('Position', [100 100 1400 600]);
    
    tl = tiledlayout(1,3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % =========================================================================
    % Density
    % =========================================================================
    
    ax1 = nexttile;
    
    imagesc(x_mm, z_mm, density_slice.');
    set(gca, 'YDir', 'reverse');
    axis image;
   
    % Remove x and y axes
    xlabel('');
    ylabel('');
    set(ax1, 'XTick', [], 'YTick', []);
    box off;
    
    % 1 cm scale bar
    hold(ax1, 'on');
    
    bar_length = 10;       % [mm]
    x0 = max(x_mm) - 13;   % position
    z0 = max(z_mm) - 2;
    
    plot(ax1, [x0 x0 + bar_length], [z0 z0], ...
        'w-', 'LineWidth', 4);
    
    text(ax1, x0 + bar_length/2, z0 - 1, '1 cm', ...
        'Color', 'w', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
    title('Density', 'FontSize', 16);

    lims = clim(ax1);
    clim(ax1, [lims(1) max_density]);
    
    cb = colorbar;
    ylabel(cb, 'Density [kg/m^3]', 'FontSize', 14);
    
    set(gca, 'FontSize', 13);
    
    % =========================================================================
    % Speed of sound
    % =========================================================================
    
    ax2 = nexttile;
    
    imagesc(x_mm, z_mm, c_slice.');
    set(gca, 'YDir', 'reverse');
    axis image;

    lims = clim(ax2);
    clim(ax2, [lims(1) max_speed]);

    % Remove x and y axes
    xlabel('');
    ylabel('');
    set(ax1, 'XTick', [], 'YTick', []);
    box off;
    
    % 1 cm scale bar
    hold(ax1, 'on');
    
    bar_length = 10;       % [mm]
    x0 = max(x_mm) - 13;   % position
    z0 = max(z_mm) - 2;
    
    plot(ax1, [x0 x0 + bar_length], [z0 z0], ...
        'w-', 'LineWidth', 4);
    
    text(ax1, x0 + bar_length/2, z0 - 1, '1 cm', ...
        'Color', 'w', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
    title('Speed of sound', 'FontSize', 16);
    
    cb = colorbar;
    ylabel(cb, 'Speed of sound [m/s]', 'FontSize', 14);
    
    set(gca, 'FontSize', 13);
    
    % =========================================================================
    % Attenuation coefficient
    % =========================================================================
    
    ax3 = nexttile;
    
    imagesc(x_mm, z_mm, alpha_slice.');
    set(gca, 'YDir', 'reverse');
    axis image;

    lims = clim(ax3);
    clim(ax3, [lims(1) max_alpha]);
    
    % Remove x and y axes
    xlabel('');
    ylabel('');
    set(ax1, 'XTick', [], 'YTick', []);
    box off;
    
    % 1 cm scale bar
    hold(ax1, 'on');
    
    bar_length = 10;       % [mm]
    x0 = max(x_mm) - 13;   % position
    z0 = max(z_mm) - 2;
    
    plot(ax1, [x0 x0 + bar_length], [z0 z0], ...
        'w-', 'LineWidth', 4);
    
    text(ax1, x0 + bar_length/2, z0 - 1, '1 cm', ...
        'Color', 'w', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');    

    title('Attenuation coefficient', 'FontSize', 16);
    
    cb = colorbar;
    ylabel(cb, '\alpha [dB/(MHz^2 cm)]', 'FontSize', 14);
    
    set(gca, 'FontSize', 13);
    
    title(tl, model_name, ...
        'FontSize', 17, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none');
    
    colormap(fig, "jet");

    % =========================================================================
    % SAVE FIGURE
    % =========================================================================
    
    [~, file_base] = fileparts(files(k).name);
    
    save_name = sprintf('AcousticProperties_%s', file_base);
    
    savefig(fig, fullfile(output_folder, [save_name '.fig']));

    exportgraphics(fig, ...
        fullfile(output_folder, [save_name '.png']), ...
        'Resolution', 300);

    clear S
end
