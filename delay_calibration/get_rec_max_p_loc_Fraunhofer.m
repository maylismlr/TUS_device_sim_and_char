%% extract_fig_peaks.m
% Loads every .fig in a folder, finds the data point with the maximum y
% value, and reports the corresponding x value, as well as saves them in .csv.
% Here to get max pressure from Fraunhofer IBMT axial pressure recordings.
%
% Author: Maylis Mueller

folder = '/data/tu_muellerm_cloud/Fraunhofer_data';
files  = dir(fullfile(folder, '*.fig'));

if isempty(files)
    error('No .fig files found in %s', folder);
end

n       = numel(files);
figName = strings(n, 1);
xAtMax  = nan(n, 1);
yMax    = nan(n, 1);
nPoints = nan(n, 1);

for k = 1:n
    fh = openfig(fullfile(folder, files(k).name), 'invisible');

    % Collect all data-bearing objects
    obj = findobj(fh, 'Type', 'line', ...
                     '-or', 'Type', 'scatter', ...
                     '-or', 'Type', 'stem');

    x = []; y = [];
    for i = 1:numel(obj)
        x = [x, obj(i).XData(:).'];
        y = [y, obj(i).YData(:).'];
    end

    figName(k) = files(k).name;

    if isempty(y)
        warning('No data objects found in %s', files(k).name);
    else
        % Ignore NaNs
        valid = ~isnan(y);
        x = x(valid);
        y = y(valid);

        [ym, idx] = max(y);
        yMax(k)    = ym;
        xAtMax(k)  = x(idx);
        nPoints(k) = numel(y);

        nTies = nnz(y == ym);
        if nTies > 1
            warning('%s: %d points tied at y = %g (first one reported)', ...
                    files(k).name, nTies, ym);
        end
    end

    close(fh);
end

T = table(figName, xAtMax, yMax, nPoints, ...
          'VariableNames', {'Figure', 'x_at_max', 'y_max', 'n_points'});

disp(T)
% writetable(T, fullfile(folder, 'fig_peaks.csv'));