%% Apply common scale to pressure field plots
% apply_common_scale_plots.m
%
% Script to load the PAC plots (or any plots) and add the same scale to
% them (up to global max p across all figures). Used to obtain some of the
% plots in Figures 3.12 and 3.13.
%
% author: Maylis Mueller

clear; clc; close all;

% =========================================================================
% FILES
% =========================================================================

fig_files = {
    % '/data/tu_muellerm/PAC/PAC_plots/MR_HR_a2862/p_fields/MR_HR_a2862_PAC_unsteered_pressured_field.fig'
    % '/data/tu_muellerm/PAC/PAC_plots/MR_HR_a2862/p_fields/MR_HR_a2862_No_PAC_unsteered_pressured_field.fig'
    % '/data/tu_muellerm/PAC/PAC_plots/MR_HR_a2862/p_fields/MR_HR_a2862_PAC_steered_pressured_field.fig'
    % '/data/tu_muellerm/PAC/PAC_plots/MR_HR_a2862/p_fields/MR_HR_a2862_No_PAC_steered_pressured_field.fig'

    '/data/tu_muellerm/PAC/PAC_plots/muCT/p_fields/muCT_PAC_unsteered_pressured_field.fig'
    '/data/tu_muellerm/PAC/PAC_plots/muCT/p_fields/muCT_No_PAC_unsteered_pressured_field.fig'
    '/data/tu_muellerm/PAC/PAC_plots/muCT/p_fields/muCT_PAC_steered_pressured_field.fig'
    '/data/tu_muellerm/PAC/PAC_plots/muCT/p_fields/muCT_No_PAC_steered_pressured_field.fig'
};

output_names = {
    % 'MR_HR_unsteered_PAC.png'
    % 'MR_HR_unsteered_No_PAC.png'
    % 'MR_HR_steered_PAC.png'
    % 'MR_HR_steered_No_PAC.png'
    'muCT_unsteered_PAC.png'
    'muCT_unsteered_No_PAC.png'
    'muCT_steered_PAC.png'
    'muCT_steered_No_PAC.png'
};

output_folder = '/data/tu_muellerm/PAC/PAC_plots/common_scale';

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

target = [0, 0, 45];   % [x y z] mm

% =========================================================================
% FIND GLOBAL MAXIMUM
% =========================================================================

global_max = 0;

for i = 1:numel(fig_files)

    fig = openfig(fig_files{i}, 'invisible');

    % Find plotted pressure image
    hImg = findobj(fig, 'Type', 'image');

    if isempty(hImg)
        error('No image object found in %s', fig_files{i});
    end

    C = hImg(1).CData;
    global_max = max(global_max, max(C(:)));

    close(fig);
end

fprintf('Common pressure scale: 0 to %.3f\n', global_max);

% =========================================================================
% APPLY COMMON SCALE AND DISPLAY
% =========================================================================

for i = 1:numel(fig_files)

    fprintf('Processing %d/%d: %s\n', ...
        i, numel(fig_files), output_names{i});

    % OPEN VISIBLE
    fig = openfig(fig_files{i}, 'visible');

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
