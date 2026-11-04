%% plot_convergence.m
% Loads and plots the convergence-study results produced by the
% convergence_study*.m scripts, for water or skull.
%
% Replaces plot_convergence_results.m and plot_convergence_results_skull.m.
%
% Medium selection (variable `medium_type` below) sets the PML result file
% and the plot style that differed between the two original scripts:
%
%                     water                       skull
%   PML file          convergence_pml_results     convergence_skull_pml_results
%   CFL file          convergence_cfl_results     convergence_skull_cfl_results
%   reference point   circle (+ x if exact zero)  x
%   CFL zero          '0' floor line              '-inf' label, no floor line
%   y upper limit     max(auto, 2 x tolerance)    fixed at 1.1
%   legend            northeast                   best
%   Goodman label     top                         bottom
%
% Study selection (variable `study` below):
%   'dx'   -> <base_dir>/dx/          convergence_dx_results.mat
%   'pml'  -> <base_dir>/PML Size/    (see table above)
%   'cfl'  -> <base_dir>/CFL number/  (see table above; full-range sweep)
%
% Each panel shows the deviation of one metric from the reference run. Runs
% that match the reference exactly have a deviation of zero; these are drawn
% on a floor below the smallest non-zero deviation.
%
% The Goodman spatial-sampling line (dx <= lambda/8) is drawn for the dx
% study only; it is a criterion on grid spacing and has no meaning on a
% PML-size or CFL axis.
%
% Figures are saved next to the results as convergence_<study>.png (water)
% or convergence_skull_<study>.png (skull).
%
% author: Maylis Mueller
% date:
% last update:

clearvars; close all;

% =========================================================================
% USER SELECTION
% =========================================================================

medium_type = 'skull';   % 'water' | 'skull'
study       = 'cfl';      % 'dx' | 'pml' | 'cfl'
base_dir    = '.';       % directory containing dx/, 'PML Size/', 'CFL number/'
save_fig    = true;

switch lower(medium_type)
    case 'water'
        is_skull = false;
        pml_file = 'convergence_pml_results.mat';
        cfl_file = 'convergence_cfl_results.mat';
    case 'skull'
        is_skull = true;
        pml_file = 'convergence_skull_pml_results.mat';
        cfl_file = 'convergence_skull_cfl_results.mat';
    otherwise
        error('medium_type must be ''water'' or ''skull''.');
end

% =========================================================================
% STUDY CONFIGURATION
% =========================================================================

switch lower(study)
    case 'dx'
        cfg = struct( ...
            'folder',  'dx', ...
            'file',    'convergence_dx_results.mat', ...
            'cand',    {{'dx'}}, ...
            'xmult',   1e3, ...        % m -> mm
            'xlabel',  'grid spacing \Deltax  [mm]', ...
            'xscale',  'log', ...
            'xdir',    'reverse', ...
            'refrule', 'min', ...      % reference = finest grid
            'goodman', true, ...
            'tag',     'dx', ...
            'title',   'grid spacing');

    case 'pml'
        cfg = struct( ...
            'folder',  'PML Size', ...
            'file',    pml_file, ...
            'cand',    {{'pml_size','pmlsize','PMLSize','pml','pml_layers'}}, ...
            'xmult',   1, ...
            'xlabel',  'PML size [grid points]', ...
            'xscale',  'linear', ...
            'xdir',    'normal', ...
            'refrule', 'max', ...      % reference = thickest PML
            'goodman', false, ...
            'tag',     'pml', ...
            'title',   'PML size');

    case 'cfl'
        cfg = struct( ...
            'folder',  'CFL number', ...
            'file',    cfl_file, ...           % full-range sweep, no merge
            'cand',    {{'cfl','CFL','cfl_number'}}, ...
            'xmult',   1, ...
            'xlabel',  'CFL number', ...
            'xscale',  'linear', ...
            'xdir',    'reverse', ...  % largest CFL on the left
            'refrule', 'min', ...      % reference = smallest time step
            'goodman', false, ...
            'tag',     'cfl', ...
            'title',   'CFL number');

    otherwise
        error('study must be ''dx'', ''pml'' or ''cfl''.');
end

% =========================================================================
% LOAD
% =========================================================================

fpath = fullfile(base_dir, cfg.folder, cfg.file);
if ~isfile(fpath)
    error('Results file not found: %s', fpath);
end

S = load(fpath);
R = S.R;
[nD, nX] = size(R);

if isfield(S, 'depth_list')
    depth_list = S.depth_list;
else
    depth_list = arrayfun(@(d) R(d,1).depth, 1:nD);
end

% --- swept parameter -----------------------------------------------------
xfield = '';
for f = cfg.cand
    if isfield(R, f{1})
        xfield = f{1};
        break;
    end
end
if isempty(xfield)
    error(['No swept-parameter field found in R.\n' ...
           'Expected one of: %s\nAvailable fields: %s'], ...
           strjoin(cfg.cand, ', '), strjoin(fieldnames(R)', ', '));
end

x = nan(nD, nX);
for d = 1:nD
    for k = 1:nX
        x(d,k) = getnum(R(d,k), xfield) * cfg.xmult;
    end
end

% --- reference case ------------------------------------------------------
if isfield(S, 'ref_idx') && ~isempty(S.ref_idx)
    ref_idx = S.ref_idx;
else
    if strcmp(cfg.refrule, 'max')
        [~, ref_idx] = max(x(1,:));
    else
        [~, ref_idx] = min(x(1,:));
    end
end

fprintf('Loaded %s\n', fpath);
fprintf('  %d depth(s), %d %s values, reference index %d (%s = %g)\n', ...
        nD, nX, xfield, ref_idx, xfield, x(1,ref_idx));

% --- axis: include all tested values, including reference ----------------
xtick    = unique(x(isfinite(x)))';
xticklab = arrayfun(@(v) num2str(v, '%.3g'), xtick, 'uni', 0);

if numel(xtick) < 2
    error('Fewer than two values available for the x-axis.');
end

if strcmp(cfg.xscale, 'log')
    xlim_err = [xtick(1)/1.08, xtick(end)*1.08];
else
    pad      = 0.06 * (xtick(end) - xtick(1));
    xlim_err = [xtick(1) - pad, xtick(end) + pad];
end

% =========================================================================
% METRIC DEFINITIONS
% =========================================================================

metrics = struct( ...
    'field',   {'pmax',            'z_from_source',         'vol_6db_iso_mm3'}, ...
    'scale',   {1e-6,              1e3,                     1}, ...
    'name',    {'peak pressure',   'axial focal position',  '-6 dB focal volume'}, ...
    'errtype', {'rel',             'abs',                   'rel'}, ...
    'errlab',  {'\epsilon_p  [%]', ...
                '\epsilon_z  [mm]', ...
                '\epsilon_V  [%]'}, ...
    'tol',     {1,                 0.5,                     1}, ...
    'tollab',  {'1 %',             '0.5 mm',                '1 %'});

nM = numel(metrics);

% =========================================================================
% DEVIATIONS
% =========================================================================

err_all = cell(1, nM);

for m = 1:nM
    E = nan(nD, nX);
    for d = 1:nD
        v = arrayfun(@(k) getnum(R(d,k), metrics(m).field), 1:nX) ...
            * metrics(m).scale;
        vref = v(ref_idx);

        switch metrics(m).errtype
            case 'rel'
                E(d,:) = abs(v - vref) ./ abs(vref) * 100;
            case 'abs'
                E(d,:) = abs(v - vref);
        end
    end
    err_all{m} = E;
end

% --- completeness report -------------------------------------------------
fprintf('\nData completeness (non-reference cases, %d per depth):\n', nX-1);
sel = true(1, nX); sel(ref_idx) = false;

for d = 1:nD
    fprintf('  %3.0f mm:', depth_list(d)*1e3);
    for m = 1:nM
        e = err_all{m}(d, sel);
        fprintf('  %s: %d finite / %d zero', ...
                metrics(m).field, sum(isfinite(e)), sum(e == 0));
    end
    fprintf('\n');
end

% =========================================================================
% FIGURE
% =========================================================================

figure('Position', [100 100 1150 350]);
cols = lines(max(nD,3));
lbl  = arrayfun(@(d) sprintf('%.0f mm', d*1e3), depth_list, 'uni', 0);

ax_err = gobjects(1, nM);

for m = 1:nM

    ax_err(m) = subplot(1, nM, m); hold on; box on;

    E  = err_all{m};
    ep = E(E > 0 & isfinite(E));

    if isempty(ep)
        yfloor = metrics(m).tol / 1e3;
    else
        yfloor = min(ep) / 3;
    end

    % exact zeros cannot sit on a log axis: clamp them to the floor
    Eplot = E;
    Eplot(E == 0) = yfloor;

    h = gobjects(1, nD);

    for d = 1:nD

        % Solid line, excluding the reference point
        Esolid = Eplot(d,:);
        Esolid(ref_idx) = NaN;

        h(d) = plot(ax_err(m), x(d,:), Esolid, 'o-', ...
            'Color', cols(d,:), 'LineWidth', 1.3, ...
            'MarkerFaceColor', 'w');

        % Point nearest to the reference
        idx_other = setdiff(1:nX, ref_idx);
        [~, ii] = min(abs(x(d,idx_other) - x(d,ref_idx)));
        idx_prev = idx_other(ii);

        % Dotted final segment to reference point
        plot(ax_err(m), ...
            [x(d,idx_prev), x(d,ref_idx)], ...
            [Eplot(d,idx_prev), Eplot(d,ref_idx)], ':', ...
            'Color', cols(d,:), ...
            'LineWidth', 1.3, ...
            'HandleVisibility', 'off');

        if is_skull
            % Reference point as cross
            plot(ax_err(m), x(d,ref_idx), Eplot(d,ref_idx), 'x', ...
                'Color', cols(d,:), ...
                'MarkerSize', 10, ...
                'LineWidth', 1.5, ...
                'HandleVisibility', 'off');

            % Only CFL: label exact zero as -infinity
            if strcmpi(study, 'cfl') && E(d,ref_idx) == 0
                text(ax_err(m), x(d,ref_idx), yfloor, '-\infty  ', ...
                    'Color', cols(d,:), ...
                    'Interpreter', 'tex', ...
                    'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'middle');
            end
        else
            % Reference point
            plot(ax_err(m), x(d,ref_idx), Eplot(d,ref_idx), 'o', ...
                'Color', cols(d,:), ...
                'MarkerFaceColor', 'w', ...
                'LineWidth', 1.3, ...
                'HandleVisibility', 'off');

            % Cross for exact zero
            if E(d,ref_idx) == 0
                plot(ax_err(m), x(d,ref_idx), yfloor, 'x', ...
                    'Color', cols(d,:), ...
                    'MarkerSize', 10, ...
                    'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    % skull CFL: log y-axis, as in convergence_study.m
    if is_skull && strcmpi(study, 'cfl')
        yscale = 'log';
    else
        yscale = 'lin';
    end

    set(ax_err(m), 'XScale', cfg.xscale, 'YScale', yscale, ...
                   'XDir', cfg.xdir, ...
                   'XTick', xtick, ...
                   'XTickLabel', xticklab, ...
                   'XLim', xlim_err, ...
                   'XMinorTick', 'off', ...
                   'XTickLabelRotation', 0);
    ax_err(m).XAxis.MinorTickValues = [];

    xlabel(ax_err(m), cfg.xlabel);
    ylabel(ax_err(m), metrics(m).errlab);
    title(ax_err(m), metrics(m).name);

    yline(ax_err(m), metrics(m).tol, 'k--', metrics(m).tollab, ...
        'LineWidth', 1.0, ...
        'LabelHorizontalAlignment', 'left', ...
        'LabelVerticalAlignment', 'bottom', ...
        'HandleVisibility', 'off');

    % zero floor line (replaced by the -inf label for skull CFL)
    if any(E(:) == 0) && ~(is_skull && strcmpi(study, 'cfl'))
        yline(ax_err(m), yfloor, 'k:', '0', ...
            'LabelHorizontalAlignment', 'right', ...
            'LabelVerticalAlignment', 'bottom', ...
            'HandleVisibility', 'off');
    end

    yl = ylim(ax_err(m));
    if is_skull && ~strcmpi(study, 'cfl')
        ylim(ax_err(m), [yfloor/2, 1.1]);
    else
        ylim(ax_err(m), [yfloor/2, max(yl(2), metrics(m).tol*2)]);
    end

    if m == 1
        if is_skull
            legend(ax_err(m), h, lbl, 'Location', 'best');
        else
            legend(ax_err(m), h, lbl, 'Location', 'northeast');
        end
    end
end

% --- Goodman criterion ---------------------------------------------------
if cfg.goodman
    if isfield(S, 'lambda') && ~isempty(S.lambda)
        lambda = S.lambda;
    else
        lambda = 1480 / 4.8e5;
    end
    dx_goodman = lambda / 8 * 1e3;   % [mm]
    fprintf('\nGoodman criterion: dx <= %.3f mm\n', dx_goodman);

    if is_skull
        goodman_valign = 'bottom';
    else
        goodman_valign = 'top';
    end

    for a = ax_err
        xline(a, dx_goodman, 'k:', 'Goodman criterion', ...
            'LineWidth', 1.3, ...
            'LabelVerticalAlignment', goodman_valign, ...
            'LabelHorizontalAlignment', 'left', ...
            'HandleVisibility', 'off');
    end
end

sgtitle(sprintf('Convergence study: %s', cfg.title));

% =========================================================================
% TABLE
% =========================================================================

fprintf('\n%-8s %-12s %-12s %-10s %-10s %-10s %-14s %-12s\n', ...
        'depth', xfield, 'pmax [MPa]', 'err_p [%]', 'z [mm]', 'dz [mm]', ...
        'V-6dB [mm3]', 'err_V [%]');
fprintf('%s\n', repmat('-', 1, 100));

for d = 1:nD
    for k = 1:nX
        fprintf('%-8.0f %-12.4g %-12.4f %-10.3f %-10.4f %-10.4f %-14.3f %-12.3f\n', ...
                depth_list(d)*1e3, x(d,k), ...
                getnum(R(d,k),'pmax')*1e-6,          err_all{1}(d,k), ...
                getnum(R(d,k),'z_from_source')*1e3,  err_all{2}(d,k), ...
                getnum(R(d,k),'vol_6db_iso_mm3'),    err_all{3}(d,k));
    end
    fprintf('\n');
end

% =========================================================================
% SAVE
% =========================================================================

if save_fig
    % skull figures get a prefix so they do not overwrite the water ones
    % in the shared study folder
    if is_skull
        figtag = ['skull_' cfg.tag];
    else
        figtag = cfg.tag;
    end
    pngfile = fullfile(base_dir, cfg.folder, ...
                       sprintf('convergence_%s.png', figtag));
    figfile = fullfile(base_dir, cfg.folder, ...
                       sprintf('convergence_%s.fig', figtag));
    exportgraphics(gcf, pngfile, 'Resolution', 200);
    savefig(gcf, figfile);
    fprintf('Saved %s\n', pngfile);
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function v = getnum(s, f)
% Scalar field access that tolerates missing or unfinished runs.
    if ~isfield(s, f) || isempty(s.(f))
        v = NaN;
    else
        v = double(s.(f)(1));
    end
end
