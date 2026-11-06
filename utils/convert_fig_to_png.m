% =========================================================================
% EXPORT ALL .FIG FILES TO PNG
% =========================================================================

folder = fullfile(cd('..'), 'simulation_plots', 'FMHW', 'sound_speed_1520');

fig_files = dir(fullfile(folder, '*.fig'));
fprintf("Number of files to convert: %d\n", length(fig_files));

for k = 1:length(fig_files)

    fprintf("Converting file %d\n", k);

    fig_name = fig_files(k).name;
    fig_path = fullfile(folder, fig_name);

    fig = openfig(fig_path, 'invisible');

    set(fig, 'PaperPositionMode', 'auto');

    [~, name, ~] = fileparts(fig_name);
    png_path = fullfile(folder, [name '.png']);

    exportgraphics(fig, png_path, 'Resolution', 300);

    close(fig);

end

disp('All .fig files exported to PNG.');