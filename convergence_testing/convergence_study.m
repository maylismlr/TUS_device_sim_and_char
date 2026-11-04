%% convergence_study.m
% Spatial (dx), PML-thickness and CFL convergence studies for the
% 256-element TRUST array model, in water or in a heterogeneous skull medium.
%
% Replaces convergence_study_dx.m, convergence_study_pml.m,
% convergence_study_cfl.m, convergence_study_skull_pml.m and
% convergence_study_skull_cfl.m. Select the case with `medium_type` and
% `study` below; every other setting follows from that choice.
%
%   study = 'dx'   sweeps the grid spacing (1.00 ... 0.20 mm) with CFL = 0.2
%                  and a 10-point PML held fixed, to demonstrate that the
%                  working resolution dx = 0.25 mm is converged. Finest
%                  grid = reference. Water only: the skull segmentation is
%                  pre-sampled at 0.25 mm and cannot follow the sweep.
%   study = 'pml'  sweeps the PML thickness (10, 20, 30 grid points) with
%                  dx = 0.25 mm and CFL = 0.2 held fixed, so that any
%                  difference between runs is attributable to boundary
%                  absorption alone. Thickest PML = reference.
%   study = 'cfl'  sweeps the CFL number (0.6 ... 0.1) with dx = 0.25 mm and
%                  a 10-point PML held fixed, so that any difference between
%                  runs is attributable to the timestep alone. The full
%                  range is run in one sweep, so no merge step is needed.
%                  Smallest CFL = reference.
%
%   medium_type = 'water'  homogeneous water, 60 x 60 x 150 mm domain,
%                          focal depths 40 / 70 / 100 mm, t_end = 160 us
%   medium_type = 'skull'  micro-CT skull in water, 56 x 56 x 70 mm domain,
%                          focal depth 40 mm, t_end = 100 us
%
% Output: results (.mat) and the study figure (.png/.fig) are written to
% the folder plot_convergence_results.m reads from:
%   <base_dir>/dx/           convergence_dx_results.mat
%   <base_dir>/PML Size/     convergence_[skull_]pml_results.mat
%   <base_dir>/CFL number/   convergence_[skull_]cfl_results.mat
%
% Metrics:
%   1) peak pressure
%   2) axial focal position
%   3) -6 dB focal volume (isosurface volume; voxel volume also saved)
%
% PMLInside = false: the layer is appended outside the physical domain, so
% changing PMLSize leaves the simulated region unaltered. This matters more
% in the skull case, because the skull is fixed to the grid and any change
% of the interior domain would shift the bone relative to the array.
%
% No calibration is applied (nominal source amplitude): this study
% characterises numerics only, so every run must sit at the same physical
% focal depth and be driven identically.
%
% NOTE: focal_metrics searches for the pressure maximum over the full
% computational domain. For skull simulations, verify that this maximum
% corresponds to the transmitted focus rather than a high-pressure region
% close to the skull or the source.
%
% author: Maylis Mueller
% date:
% last update:

clearvars;
close all;

% =========================================================================
% CASE SELECTION
% =========================================================================

medium_type = 'skull';   % 'water' | 'skull'
study       = 'pml';      % 'dx' | 'pml' | 'cfl'
base_dir    = '.';       % directory containing dx/, 'PML Size/', 'CFL number/'

medium_type = lower(medium_type);
study       = lower(study);
key         = [medium_type '_' study];

if strcmp(key, 'skull_dx')
    error(['No dx study for the skull: the segmentation is pre-sampled at ' ...
           '0.25 mm and is not resampled to other grid spacings.']);
end
if ~any(strcmp(key, {'water_dx', 'water_pml', 'water_cfl', 'skull_pml', 'skull_cfl'}))
    error('medium_type must be ''water'' or ''skull''; study must be ''dx'', ''pml'' or ''cfl''.');
end

addpath '/data/hu_muellerm/MATLAB Add-Ons/Collections/k-Wave/k-Wave'
addpath '/data/tu_muellerm/utils';

% output folder, matching plot_convergence_results.m
switch study
    case 'dx',  folder = 'dx';
    case 'pml', folder = 'PML Size';
    case 'cfl', folder = 'CFL number';
end
outdir = fullfile(base_dir, folder);
if ~isfolder(outdir)
    mkdir(outdir);
end

% Results are written to a file of their own per case: overwriting another
% sweep would silently mix studies with incompatible result structures.
switch key
    case 'water_dx',  fname = 'convergence_dx_results.mat';
    case 'water_pml', fname = 'convergence_pml_results.mat';
    case 'water_cfl', fname = 'convergence_cfl_results.mat';
    case 'skull_pml', fname = 'convergence_skull_pml_results.mat';
    case 'skull_cfl', fname = 'convergence_skull_cfl_results.mat';
end
outfile = fullfile(outdir, fname);

% =========================================================================
% USER-DEFINED LITERALS
% =========================================================================

% medium (water background)
c0   = 1480;           % sound speed [m/s]
rho0 = 1000;           % density [kg/m^3]

% source
source_f0     = 480e3;   % centre frequency [Hz]
source_cycles = 30;      % 62.5 us burst at 480 kHz
num_bursts    = 1;
PRF_period    = 15e-6;   % [s]

element_width  = 2.58e-3;        % width [m]
element_height = element_width;  % height [m]
element_pitch  = 3.08e-3;        % pitch [m]

row_counts = [6 10 12 14 16 16 18 18 18 18 18 18 16 16 14 12 10 6];

% wavelength in water; water gives the smallest wavelength in this model
lambda = c0 / source_f0;   % 3.083 mm

% working spatial discretisation (overridden per run in the dx study)
dx = 0.25e-3;          % [m]
dy = dx;
dz = dx;

switch medium_type
    case 'water'
        source_amp = 0.1e6;    % [Pa] arbitrary but IDENTICAL across all runs

        % fixed physical domain (identical for every run)
        Lx = 60e-3;
        Ly = 60e-3;
        Lz = 150e-3;

        t_end      = 160e-6;               % simulation time [s]
        depth_list = [40, 70, 100] * 1e-3; % focal depths [m]

    case 'skull'
        % skull segmentation (micro-CT, pre-sampled at 0.25 mm)
        path_medium = '/data/tu_muellerm_cloud/skulls_FUN/skulls_params/final/muCT.nrrd_0.25mm_params.mat';
        tmp = load(path_medium);
        seg = tmp.seg;
        clear tmp;

        source_amp = 1.5e6;    % [Pa] arbitrary but IDENTICAL across all runs

        % physical computational domain
        Lx = 56e-3;
        Ly = 56e-3;
        Lz = 70e-3;

        t_end      = 100e-6;   % simulation time [s]
        depth_list = 40e-3;    % focal depth [m]
end

% =========================================================================
% SWEEP DEFINITION
% =========================================================================

switch study
    case 'dx'
        cfl        = 0.2;            % fixed CFL number
        pml_size   = 10;             % fixed PML thickness [grid points]
        dx_list    = [1.00 0.75 0.50 0.40 0.30 0.25 0.20] * 1e-3;   % [m]
        sweep_list = dx_list;
        list_name  = 'dx_list';

        % finest grid = reference
        ref_idx = numel(dx_list);

    case 'pml'
        cfl        = 0.2;            % fixed CFL number
        pml_list   = [10 20 30];     % PML thickness [grid points]
        sweep_list = pml_list;
        list_name  = 'pml_list';

        % thickest absorbing layer = reference
        [~, ref_idx] = max(pml_list);

    case 'cfl'
        pml_size   = 10;             % fixed PML thickness [grid points]
        cfl_list   = [0.6 0.5 0.4 0.3 0.2 0.1];
        sweep_list = cfl_list;
        list_name  = 'cfl_list';

        % finest temporal resolution = reference
        [~, ref_idx] = min(cfl_list);
end

% Field layout of R and the variables written to outfile are kept as in the
% original per-case scripts, so plot_convergence_results.m reads them
% unchanged.
R_common = {'dx',[], 'ppw',[], 'Nx',[], 'Ny',[], 'Nz',[], 'dt',[], 'Nt',[], ...
            'pmax',[], 'z_from_source',[], 'x_pos',[], 'y_pos',[], ...
            'vol_6db_voxel_mm3',[], 'vol_6db_iso_mm3',[], ...
            'runtime',[], 'done',false};

switch key
    case 'water_dx'
        xfield     = 'dx';
        R_template = struct('depth',[], R_common{:});
        save_vars  = {'R', 'dx_list', 'depth_list', 'ref_idx', 'lambda'};
    case 'water_pml'
        xfield     = 'pml';
        R_template = struct('depth',[], 'pml',[], R_common{:});
        save_vars  = {'R', 'pml_list', 'depth_list', 'ref_idx', 'lambda', ...
                      'dx', 'c0', 'source_f0'};
    case 'skull_pml'
        xfield     = 'pml_size';
        R_template = struct('depth',[], 'pml_size',[], 'cfl',[], R_common{:});
        save_vars  = {'R', 'pml_list', 'depth_list', 'ref_idx', 'lambda', ...
                      'dx', 'cfl', 'c0', 'source_f0'};
    otherwise   % water_cfl, skull_cfl
        xfield     = 'cfl';
        R_template = struct('depth',[], 'cfl',[], R_common{:});
        save_vars  = {'R', 'cfl_list', 'depth_list', 'ref_idx', 'lambda', ...
                      'dx', 'c0', 'source_f0'};
end

% =========================================================================
% GRID
% =========================================================================

% In the dx study the grid is rebuilt inside the loop for every spacing.
if ~strcmp(study, 'dx')
    Nx = roundEven(Lx / dx);
    Ny = roundEven(Ly / dy);
    Nz = roundEven(Lz / dz);

    switch medium_type
        case 'water'
            % Only dt (CFL study) or the externally appended PML (PML
            % study) varies across the sweep.
            fprintf('Fixed grid: %d x %d x %d\n', Nx, Ny, Nz);
            fprintf('dx = %.3f mm, PPW = %.2f\n', dx*1e3, lambda/dx);
        case 'skull'
            fprintf('Grid: %d x %d x %d\n', Nx, Ny, Nz);
            fprintf('dx = %.3f mm\n', dx*1e3);
            fprintf('PPW in water = %.2f\n', lambda/dx);
    end
end

% =========================================================================
% MEDIUM
% =========================================================================

switch medium_type
    case 'water'
        medium.sound_speed = c0;
        medium.density     = rho0;

    case 'skull'
        % The segmentation is exported on its own voxel grid, so it is
        % checked against the simulation grid rather than resampled: a
        % silent mismatch here would displace the skull relative to the
        % array.
        fprintf('Loaded segmentation: %d x %d x %d voxels\n', size(seg,1), size(seg,2), size(seg,3));

        % lateral extent must match exactly
        if size(seg,1) ~= Nx || size(seg,2) ~= Ny
            error(['Segmentation x-y dimensions do not match the simulation grid.\n' ...
                   'seg  = %d x %d\ngrid = %d x %d'], size(seg,1), size(seg,2), Nx, Ny);
        end

        % axial extent may exceed the domain, but not fall short of it
        if size(seg,3) < Nz
            error('Segmentation contains only %d z-slices, but the simulation requires %d.', ...
                  size(seg,3), Nz);
        end

        % crop to the simulation domain and enforce a binary mask
        seg = seg(:,:,1:Nz);
        seg = seg > 0.5;

        fprintf('Segmentation used: %d x %d x %d voxels\n', size(seg,1), size(seg,2), size(seg,3));

        medium.sound_speed = zeros(Nx, Ny, Nz);
        medium.density     = zeros(Nx, Ny, Nz);
        medium.alpha_coeff = zeros(Nx, Ny, Nz);

        % water
        medium.sound_speed(seg == 0) = 1480;
        medium.density(seg == 0)     = 1000;
        medium.alpha_coeff(seg == 0) = 2.17e-3;    % dB/(MHz^2 cm)

        % skull
        medium.sound_speed(seg == 1) = 3100;
        medium.density(seg == 1)     = 2000;
        medium.alpha_coeff(seg == 1) = 2.7;        % dB/(MHz^2 cm)

        % y = 2 makes the k-Wave dispersion term vanish
        medium.alpha_power = 2;

        fprintf('Maximum sound speed = %.0f m/s\n', max(medium.sound_speed(:)));
end

% =========================================================================
% SWEEP
% =========================================================================

nD = numel(depth_list);
nX = numel(sweep_list);

% Resume if a partial run exists, but only after confirming that the stored
% sweep matches the one requested here. A sweep or depth list that has been
% edited since the last run invalidates the saved results entirely.
resume_ok = false;

if isfile(outfile)
    S = load(outfile);
    if isfield(S,'R') && isfield(S,list_name) && isfield(S,'depth_list') && ...
       isequal(S.(list_name), sweep_list) && isequal(S.depth_list, depth_list) && ...
       isequal(size(S.R), [nD, nX])
        R = S.R;
        resume_ok = true;
        fprintf('Resuming from %s\n', outfile);
    else
        warning(['Existing %s does not correspond to the current %s sweep. ' ...
                 'Starting a new result structure.'], outfile, upper(study));
    end
end

if ~resume_ok
    R = repmat(R_template, nD, nX);
end

for d = 1:nD
    depth = depth_list(d);

    for k = 1:nX

        % swept value for this run
        switch study
            case 'dx'
                dx   = dx_list(k);
                dy   = dx;
                dz   = dx;
                xval = dx;
            case 'pml'
                pml_size = pml_list(k);
                xval     = pml_size;
            case 'cfl'
                cfl  = cfl_list(k);
                xval = cfl;
        end

        % console labels, as printed by the original scripts
        skip_note = 'already completed';
        switch key
            case 'water_dx'
                skip_tag  = sprintf('dx %.3f mm', dx*1e3);
                skip_note = 'already done with volume';
            case 'water_pml'
                skip_tag = sprintf('PML size %.3f', pml_size);
                run_tag  = sprintf('PML %.3f', pml_size);
            case 'skull_pml'
                skip_tag = sprintf('PML %d points', pml_size);
                run_tag  = sprintf('PML %2d points | CFL %.2f', pml_size, cfl);
            otherwise
                skip_tag = sprintf('CFL %.3f', cfl);
                run_tag  = skip_tag;
        end

        % skip completed cases
        has_volume = ~isempty(R(d,k).vol_6db_iso_mm3) && isfinite(R(d,k).vol_6db_iso_mm3);
        same_x     = ~isempty(R(d,k).(xfield)) && abs(R(d,k).(xfield) - xval) < 1e-12;

        if R(d,k).done && has_volume && same_x
            fprintf('skip: depth %.0f mm, %s (%s)\n', depth*1e3, skip_tag, skip_note);
            continue;
        end

        tic;

        % grid (dx study only: rebuilt for every spacing)
        if strcmp(study, 'dx')
            Nx = roundEven(Lx / dx);
            Ny = roundEven(Ly / dy);
            Nz = roundEven(Lz / dz);
        end

        % grid and timestep
        % dt is set from the maximum sound speed in the medium (c0 in water,
        % the skull speed in the heterogeneous case).
        kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
        kgrid.makeTime(medium.sound_speed, cfl, t_end);

        % focus, defined relative to the array face
        focal_point_grid = [0, 0, kgrid.z_vec(1) + depth];

        % source
        % Delays are computed for the water sound speed: the beamformer is
        % not given the skull, so the geometric focus is the commanded one.
        source = create_source(focal_point_grid, row_counts, element_height, ...
                               element_width, element_pitch, source_amp, ...
                               1/kgrid.dt, source_f0, c0, kgrid, ...
                               'NumCycles', source_cycles);

        % sensor
        clear sensor;
        sensor.mask   = true(Nx, Ny, Nz);
        sensor.record = {'p_max'};

        % run
        sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, ...
            'PMLSize',     pml_size, ...
            'PMLInside',   false, ...
            'PlotPML',     false, ...
            'DisplayMask', 'off', ...
            'PlotSim',     false);

        P = reshape(sensor_data.p_max, Nx, Ny, Nz);

        % peak pressure, focal position and -6 dB focal volume
        % M.idx supplies the focal peak explicitly, so the volume extraction
        % is not misled by supra-threshold pressure inside the bone and is
        % independent of resolution-dependent absolute voxel indices.
        % (The dx study keeps its original focal_vol_dB plot flag = true.)
        M = focal_metrics(P, kgrid);
        V = focal_vol_dB(P, dx, dy, dz, strcmp(study, 'dx'), M.idx);

        % store
        R(d,k).depth = depth;
        if strcmp(study, 'pml')
            R(d,k).(xfield) = pml_size;
        end
        if isfield(R, 'cfl')
            R(d,k).cfl = cfl;
        end
        R(d,k).dx    = dx;
        R(d,k).ppw   = lambda / dx;
        R(d,k).Nx    = Nx;
        R(d,k).Ny    = Ny;
        R(d,k).Nz    = Nz;
        R(d,k).dt    = kgrid.dt;
        R(d,k).Nt    = kgrid.Nt;
        R(d,k).pmax  = M.pmax;
        R(d,k).z_from_source = M.pos(3) - kgrid.z_vec(1);
        R(d,k).x_pos = M.pos(1);
        R(d,k).y_pos = M.pos(2);
        R(d,k).vol_6db_voxel_mm3 = V.V_focal_mm3;
        R(d,k).vol_6db_iso_mm3   = V.V_isosurf_mm3;
        if strcmp(study, 'dx')
            R(d,k).vol_used_peak = V.p_peak;
        end
        R(d,k).runtime = toc;
        R(d,k).done    = true;

        if strcmp(study, 'dx')
            fprintf(['depth %3.0f mm | dx %.3f mm (PPW %5.2f) | %4dx%4dx%4d | ' ...
                     'pmax %.4f MPa | z %.3f mm | V-6dB %.2f mm^3 | %.0f s\n'], ...
                    depth*1e3, dx*1e3, lambda/dx, Nx, Ny, Nz, ...
                    M.pmax*1e-6, R(d,k).z_from_source*1e3, ...
                    R(d,k).vol_6db_iso_mm3, R(d,k).runtime);
        else
            fprintf(['depth %3.0f mm | %s | dt %6.2f ns | Nt %5d | ' ...
                     'pmax %.4f MPa | z %.3f mm | V-6dB %.2f mm^3 | %.0f s\n'], ...
                    depth*1e3, run_tag, kgrid.dt*1e9, kgrid.Nt, M.pmax*1e-6, ...
                    R(d,k).z_from_source*1e3, R(d,k).vol_6db_iso_mm3, R(d,k).runtime);
        end

        % save after every completed run
        save(outfile, save_vars{:});

        clear P sensor_data source kgrid M V;
    end
end

% =========================================================================
% ERROR RELATIVE TO THE REFERENCE RUN
% =========================================================================

err_p_all = nan(nD, nX);
err_z_all = nan(nD, nX);
err_V_all = nan(nD, nX);

for d = 1:nD
    ref = R(d, ref_idx);
    for k = 1:nX
        err_p_all(d,k) = abs(R(d,k).pmax - ref.pmax) / abs(ref.pmax) * 100;
        err_z_all(d,k) = abs(R(d,k).z_from_source - ref.z_from_source) * 1e3;
        err_V_all(d,k) = abs(R(d,k).vol_6db_iso_mm3 - ref.vol_6db_iso_mm3) / ...
                         abs(ref.vol_6db_iso_mm3) * 100;
    end
end

fprintf('\n');
switch key
    case 'water_dx'
        % (no reference line in the original dx script)
    case 'water_pml', fprintf('Reference PML size = %.3f\n\n', pml_list(ref_idx));
    case 'skull_pml', fprintf('Reference PML = %d grid points\n\n', pml_list(ref_idx));
    otherwise,        fprintf('Reference CFL = %.3f\n\n', cfl_list(ref_idx));
end

switch key
    case 'water_dx'
        fprintf('%-8s %-8s %-7s %-12s %-10s %-10s %-14s %-12s\n', ...
                'depth', 'dx', 'PPW', 'pmax [MPa]', 'err_p [%]', 'dz [mm]', ...
                'V-6dB [mm3]', 'err_V [%]');
        fprintf('%s\n', repmat('-', 1, 98));

        for d = 1:nD
            for k = 1:nX
                fprintf('%-8.0f %-8.3f %-7.2f %-12.4f %-10.3f %-10.4f %-14.3f %-12.3f\n', ...
                        R(d,k).depth*1e3, R(d,k).dx*1e3, R(d,k).ppw, ...
                        R(d,k).pmax*1e-6, err_p_all(d,k), err_z_all(d,k), ...
                        R(d,k).vol_6db_iso_mm3, err_V_all(d,k));
            end
            fprintf('\n');
        end

    case 'skull_pml'
        fprintf('%-8s %-10s %-10s %-12s %-10s %-10s %-14s %-12s\n', ...
                'depth', 'PML', 'dt [ns]', 'pmax [MPa]', 'err_p [%]', 'dz [mm]', ...
                'V-6dB [mm3]', 'err_V [%]');
        fprintf('%s\n', repmat('-', 1, 100));

        for d = 1:nD
            for k = 1:nX
                fprintf('%-8.0f %-10d %-10.3f %-12.4f %-10.3f %-10.4f %-14.3f %-12.3f\n', ...
                        R(d,k).depth*1e3, R(d,k).pml_size, R(d,k).dt*1e9, ...
                        R(d,k).pmax*1e-6, err_p_all(d,k), err_z_all(d,k), ...
                        R(d,k).vol_6db_iso_mm3, err_V_all(d,k));
            end
            fprintf('\n');
        end

    otherwise   % water_pml, water_cfl, skull_cfl
        fprintf('%-8s %-8s %-10s %-8s %-12s %-10s %-10s %-14s %-12s\n', ...
                'depth', upper(study), 'dt [ns]', 'Nt', 'pmax [MPa]', 'err_p [%]', ...
                'dz [mm]', 'V-6dB [mm3]', 'err_V [%]');
        fprintf('%s\n', repmat('-', 1, 105));

        for d = 1:nD
            for k = 1:nX
                fprintf('%-8.0f %-8.3f %-10.3f %-8d %-12.4f %-10.3f %-10.4f %-14.3f %-12.3f\n', ...
                        R(d,k).depth*1e3, R(d,k).(xfield), R(d,k).dt*1e9, R(d,k).Nt, ...
                        R(d,k).pmax*1e-6, err_p_all(d,k), err_z_all(d,k), ...
                        R(d,k).vol_6db_iso_mm3, err_V_all(d,k));
            end
            fprintf('\n');
        end
end

% =========================================================================
% FIGURE
% =========================================================================

lbl = arrayfun(@(d) sprintf('%.0f mm', d*1e3), depth_list, 'uni', 0);

% panel definitions
E_all       = {err_p_all, err_z_all, err_V_all};
plot_titles = {'peak pressure', 'axial focal position', '-6 dB focal volume'};
ylabs       = {'|\Delta p_{max}| / p_{max}^{ref} [%]', '|\Delta z| [mm]', ...
               '|\Delta V_{-6dB}| / V_{-6dB}^{ref} [%]'};
tols        = [1, 0.5, 1];
tollabs     = {'1 %', '0.5 mm', '1 %'};

switch study
    case 'dx',  xlab = 'points per wavelength';  sg = 'grid spacing';
    case 'pml', xlab = 'PML size [grid points]'; sg = 'PML size';
    case 'cfl', xlab = 'CFL number';             sg = 'CFL number';
end

switch key

    % ---------------------------------------------------------------------
    % water, dx: error vs points per wavelength, Goodman line
    % ---------------------------------------------------------------------
    case 'water_dx'
        % Goodman spatial sampling criterion
        dx_goodman  = lambda / 8;
        ppw_goodman = lambda / dx_goodman;

        fprintf('Goodman criterion: dx <= %.3f mm (%.2f PPW)\n', ...
                dx_goodman*1e3, ppw_goodman);

        figure('Position', [100 100 1100 350]);
        cols = lines(nD);

        for m = 1:3
            subplot(1,3,m); hold on; box on;

            % the reference (zero error) drops out of the log axis
            for d = 1:nD
                plot([R(d,:).ppw], E_all{m}(d,:), 'o-', ...
                     'Color', cols(d,:), 'LineWidth', 1.3);
            end

            set(gca, 'XScale', 'log', 'YScale', 'log');
            yline(tols(m), 'k--');
            xline(ppw_goodman, 'k:', 'Goodman criterion', ...
                  'LineWidth', 1.3, ...
                  'LabelVerticalAlignment', 'bottom', ...
                  'LabelHorizontalAlignment', 'left');

            xlabel(xlab);
            ylabel(ylabs{m});
            title(plot_titles{m});
            if m == 1
                legend(lbl, 'Location', 'southwest');
            end
        end

    % ---------------------------------------------------------------------
    % water, PML: zeros clamped to a fixed floor, fixed y-limits
    % ---------------------------------------------------------------------
    case 'water_pml'
        figure('Position', [100 100 1100 350]);
        cols = lines(nD);

        % Lowest value that can be displayed on the logarithmic axis.
        % Exact zero errors are plotted here for visualization only.
        floors = [1e-3, 1e-4, 1e-3];                 % [%], [mm], [%]
        ylims  = {[1e-3 100], [1e-4 10], [1e-3 100]};

        for m = 1:3
            subplot(1,3,m); hold on; box on;

            for d = 1:nD
                plot([R(d,:).pml], max(E_all{m}(d,:), floors(m)), 'o-', ...
                     'Color', cols(d,:), 'LineWidth', 1.3);
            end

            % only three thicknesses are sampled, so the x-axis stays
            % linear and is ticked at the sampled values
            set(gca, 'YScale', 'log', 'XTick', pml_list);
            ylim(ylims{m});
            yline(tols(m), 'k--');

            xlabel(xlab);
            ylabel(ylabs{m});
            title(plot_titles{m});
            if m == 1
                legend(lbl, 'Location', 'best');
            end
        end

    % ---------------------------------------------------------------------
    % water, CFL: reference omitted, log-log, reversed x-axis
    % ---------------------------------------------------------------------
    case 'water_cfl'
        figure('Position', [100 100 1100 350]);
        cols = lines(nD);

        for m = 1:3
            subplot(1,3,m); hold on; box on;

            for d = 1:nD
                e = E_all{m}(d,:);

                % Reference has zero error and cannot be shown on log y-axis
                e(ref_idx) = NaN;

                plot([R(d,:).cfl], e, 'o-', 'Color', cols(d,:), 'LineWidth', 1.3);
            end

            % axis reversed so that resolution increases to the right, as in
            % the spatial convergence figure
            set(gca, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse');
            yline(tols(m), 'k--');

            xlabel(xlab);
            ylabel(ylabs{m});
            title(plot_titles{m});
            if m == 1
                legend(lbl, 'Location', 'best');
            end
        end

    % ---------------------------------------------------------------------
    % skull, PML or CFL: zeros at a data-driven floor, marked with an x
    % ---------------------------------------------------------------------
    otherwise
        figure('Position', [100 100 1150 350]);
        cols = lines(max(nD,3));

        % x-axis: all sampled values are ticked, reference included; the
        % axis is linear, and reversed for CFL so that temporal resolution
        % increases to the right.
        if strcmp(study, 'cfl')
            xdir = 'reverse';
        else
            xdir = 'normal';
        end
        xtick    = sort(unique(sweep_list));
        pad      = 0.06 * (max(xtick) - min(xtick));
        xlim_err = [min(xtick) - pad, max(xtick) + pad];

        ax = gobjects(1,3);

        for m = 1:3
            ax(m) = subplot(1,3,m);
            hold(ax(m), 'on');
            box(ax(m), 'on');

            E = E_all{m};

            % The reference has exactly zero error. Zero cannot be shown on
            % a logarithmic axis, so it is displayed at a visual floor and
            % marked with an x.
            ep = E(E > 0 & isfinite(E));
            if isempty(ep)
                yfloor = tols(m) / 1e3;
            else
                yfloor = min(ep) / 3;
            end

            % one line per focal depth
            h = gobjects(1,nD);

            for d = 1:nD
                Eplot = E(d,:);
                Eplot(Eplot == 0) = yfloor;

                h(d) = plot(ax(m), sweep_list, Eplot, 'o-', ...
                            'Color', cols(d,:), 'LineWidth', 1.3, 'MarkerFaceColor', 'w');

                z = find(E(d,:) == 0);
                if ~isempty(z)
                    plot(ax(m), sweep_list(z), yfloor * ones(size(z)), 'x', ...
                         'Color', cols(d,:), 'MarkerSize', 10, 'LineWidth', 1.5, ...
                         'HandleVisibility', 'off');
                end
            end

            % axes
            set(ax(m), 'XScale','linear', 'YScale','log', 'XDir', xdir, ...
                       'XTick', xtick, 'XLim', xlim_err, 'XMinorTick','off', ...
                       'XTickLabelRotation', 0);

            xlabel(ax(m), xlab);
            ylabel(ax(m), ylabs{m});
            title(ax(m), plot_titles{m});

            % convergence criterion
            yline(ax(m), tols(m), 'k--', tollabs{m}, 'LineWidth', 1.0, ...
                  'LabelHorizontalAlignment', 'left', ...
                  'LabelVerticalAlignment', 'bottom', ...
                  'HandleVisibility', 'off');

            % y limits, bounded below by the visual floor rather than by an
            % artificial zero line, and above by at least twice the tolerance
            yl = ylim(ax(m));
            ylim(ax(m), [yfloor/2, max(yl(2), tols(m)*2)]);

            if m == 1
                legend(ax(m), h, lbl, 'Location', 'best');
            end
        end

        sgtitle(['Convergence study: ' sg]);
end

% -------------------------------------------------------------------------
% save (next to the results, in the study folder)
% -------------------------------------------------------------------------
% convergence_dx_study, convergence_pml_study, convergence_cfl_study,
% convergence_skull_pml_study, convergence_skull_cfl_study
figbase = strrep(outfile, '_results.mat', '_study');
exportgraphics(gcf, [figbase '.png'], 'Resolution', 200);
savefig(gcf, [figbase '.fig']);
fprintf('Saved %s\n%s\n', outfile, [figbase '.png']);
