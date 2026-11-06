function volume = focal_vol_dB(p_max_3D, dx, dy, dz, do_plot, focus_idx, varargin)
% FOCAL_VOL_DB  Focal volume of a 3D peak-pressure field at a configurable dB level.
%
% POSITIONAL INPUTS:
%   p_max_3D  : 3D array of peak pressure, dimensions ordered (x, y, z)
%   dx,dy,dz  : voxel spacing (m) along dimensions 1, 2, 3
%   do_plot   : (optional) true/false to show isosurface plot. Default false.
%   focus_idx : (optional) [ix iy iz] index of the focal peak. When supplied,
%               this voxel defines the focal peak and connected component.
%               This is the recommended mode for convergence testing because
%               it avoids any resolution-dependent hard-coded z index.
%
% NAME-VALUE INPUTS:
%   'LevelDB'      : threshold level in dB below the focal peak. Default -6.
%                      -6 dB -> 50.0 % of peak pressure (25 % of peak intensity)
%                               = FWHM of the pressure profile
%                      -3 dB -> 70.7 % of peak pressure (50 % of peak intensity)
%                               = FWHM of the intensity profile
%                    NOTE: a dB level is unambiguous. It is the same number
%                    whether computed as 20*log10(p/p0) or 10*log10(I/I0).
%                    The -3 dB and -6 dB contours are therefore NOT the same
%                    surface, and differ in linear width by ~sqrt(2).
%   'DBConvention' : 'half' (default) | 'strict'
%                      'half'   -> canonical rounded fractions, p_frac = 2^(L/6):
%                                  -6 dB = 0.5, -3 dB = 1/sqrt(2), -12 dB = 0.25.
%                                  This is what the literature means by "-6 dB".
%                      'strict' -> literal p_frac = 10^(L/20) (-6 dB = 0.50119).
%                    The two differ by <0.5 % in volume; 'half' is the default
%                    for consistency with ITRUSST / IEC reporting practice.
%   'Verbose'      : print threshold and results. Default false (quiet in loops).
%
% OUTPUT (struct):
%   volume.N_voxels       : number of voxels in the thresholded focal region
%   volume.V_focal_mm3    : voxel-counting volume (mm^3)
%   volume.V_isosurf_mm3  : marching-cubes isosurface volume (mm^3)
%   volume.p_peak         : pressure value defining the threshold
%   volume.p_thresh       : absolute threshold pressure
%   volume.p_frac         : threshold as a fraction of p_peak
%   volume.level_dB       : threshold level actually applied
%   volume.focus_idx      : [ix iy iz] used as the focal peak
%   volume.truncated      : true if the focal component touches a domain face,
%                           in which case V_isosurf_mm3 is not trustworthy
%
%   author: Maylis Mueller

    % ---------------- Option parsing ----------------
    if nargin < 5 || isempty(do_plot)
        do_plot = false;
    end
    if nargin < 6
        focus_idx = [];
    end

    opt.LevelDB      = -6;
    opt.DBConvention = 'half';
    opt.Verbose      = false;

    if mod(numel(varargin), 2) ~= 0
        error('Name-value arguments must come in pairs.');
    end
    for a = 1:2:numel(varargin)
        name = varargin{a};
        if ~isfield(opt, name)
            error('Unknown option "%s".', name);
        end
        opt.(name) = varargin{a+1};
    end

    if opt.LevelDB > 0
        error('LevelDB must be <= 0 (it is a level below the focal peak).');
    end

    switch lower(opt.DBConvention)
        case 'half'
            p_frac = 2^(opt.LevelDB / 6);       % -6 -> 0.5, -3 -> 1/sqrt(2)
        case 'strict'
            p_frac = 10^(opt.LevelDB / 20);     % -6 -> 0.50119
        otherwise
            error('DBConvention must be ''half'' or ''strict''.');
    end

    % ---------------- Define focal peak ----------------
    % NEEDS ACTUAL RECORDED PEAK IN SIM
    if ~isempty(focus_idx)
        ix0 = focus_idx(1);
        iy0 = focus_idx(2);
        iz0 = focus_idx(3);
        p_max = p_max_3D(ix0, iy0, iz0);
    else
        % Backwards-compatible behaviour for previous skull analyses.
        % NOTE: for convergence testing, always pass focus_idx.
        if size(p_max_3D, 3) < 80
            error(['p_max_3D has fewer than 80 z samples. Pass focus_idx to ', ...
                   'focal_vol_dB instead of using the legacy z>=80 search.']);
        end
        sub_vol = p_max_3D(:,:,80:end);
        p_max = max(sub_vol(:));
        [ix0, iy0, iz0_sub] = ind2sub(size(sub_vol), find(sub_vol == p_max, 1));
        iz0 = iz0_sub + 79;
    end

    if ~isfinite(p_max) || p_max <= 0
        error('Focal peak pressure is %g. Check p_max_3D and focus_idx.', p_max);
    end

    volume.p_peak    = p_max;
    volume.p_frac    = p_frac;
    volume.level_dB  = opt.LevelDB;
    volume.focus_idx = [ix0 iy0 iz0];

    % ---------------- Threshold and connected component ----------------
    p_thresh = p_max * p_frac;
    volume.p_thresh = p_thresh;

    mask = p_max_3D >= p_thresh;

    CC = bwconncomp(mask, 26);
    focus_idx_linear = sub2ind(size(p_max_3D), ix0, iy0, iz0);

    focus_label = 0;
    for k = 1:CC.NumObjects
        if any(CC.PixelIdxList{k} == focus_idx_linear)
            focus_label = k;
            break;
        end
    end

    if focus_label == 0
        warning('Focus voxel not found in any connected component. Check p_max_3D indexing.');
        volume.N_voxels      = NaN;
        volume.V_focal_mm3   = NaN;
        volume.V_isosurf_mm3 = NaN;
        volume.truncated     = true;
        return;
    end

    focus_pixels = CC.PixelIdxList{focus_label};

    % ---------------- Voxel-counting volume ----------------
    voxel_volume       = dx * dy * dz;
    volume.N_voxels    = numel(focus_pixels);
    volume.V_focal_mm3 = volume.N_voxels * voxel_volume * 1e9;

    % ---------------- Isosurface (marching cubes) volume ----------------
    [ii, jj, kk] = ind2sub(size(mask), focus_pixels);

    volume.truncated = min(ii) <= 1 || max(ii) >= size(mask,1) || ...
                       min(jj) <= 1 || max(jj) >= size(mask,2) || ...
                       min(kk) <= 1 || max(kk) >= size(mask,3);
    if volume.truncated
        warning(['Focal component touches a domain face. The isosurface is ', ...
                 'closed artificially and V_isosurf_mm3 is not trustworthy.']);
    end

    % Marching cubes runs on the full field, wrapped in a single sub-threshold
    % layer so the surface is always closed (isosurface does not cap at array
    % boundaries, and an open mesh makes the divergence-theorem volume
    % meaningless). Allocating the padded array and copying into it once avoids
    % holding a second full-size working copy for the masking step below.
    p_pad = zeros(size(p_max_3D) + 2, 'like', p_max_3D);
    p_pad(2:end-1, 2:end-1, 2:end-1) = p_max_3D;

    % Suppress supra-threshold voxels belonging to OTHER components so their
    % closed surfaces do not contribute to the volume sum. The true field
    % values are left everywhere else: marching cubes places each vertex at
    % t = (iso - a)/(b - a) along the edge between an inside voxel a and an
    % outside voxel b, so blanking b would pull every boundary vertex onto the
    % last inside voxel and bias the volume low by ~half a voxel over the whole
    % surface, by an amount that scales with dx. Under 26-connectivity, voxels
    % of other components are never adjacent to the focal component, so zeroing
    % them cannot perturb any edge crossing on the focal boundary.
    other_labels = setdiff(1:CC.NumObjects, focus_label);
    if ~isempty(other_labels)
        other_pixels = vertcat(CC.PixelIdxList{other_labels});
        [oi, oj, ok] = ind2sub(size(mask), other_pixels);
        p_pad(sub2ind(size(p_pad), oi + 1, oj + 1, ok + 1)) = 0;
    end

    fv_focal = isosurface(p_pad, p_thresh);
    clear p_pad

    if isempty(fv_focal.vertices)
        volume.V_isosurf_mm3 = NaN;
    else
        % With isosurface(V,...), MATLAB vertex column 1 corresponds to array
        % dimension 2, column 2 to dimension 1, and column 3 to dimension 3.
        % Translation does not affect enclosed volume, so only the spacing
        % scale factors are needed here.
        vertices_m = fv_focal.vertices;
        vertices_m(:,1) = vertices_m(:,1) * dy;
        vertices_m(:,2) = vertices_m(:,2) * dx;
        vertices_m(:,3) = vertices_m(:,3) * dz;

        volume.V_isosurf_mm3 = isosurfaceVolume(vertices_m, fv_focal.faces) * 1e9;
    end

    % ---------------- Reporting ----------------
    if opt.Verbose
        fprintf('Threshold: %.1f dB (%s) -> %.4f x p_peak = %.4f MPa\n', ...
                opt.LevelDB, lower(opt.DBConvention), p_frac, p_thresh * 1e-6);
        fprintf('Focal peak: %.4f MPa at index [%d %d %d]\n', ...
                p_max * 1e-6, ix0, iy0, iz0);
        fprintf('Focal volume (voxel count): %.4f mm^3  (%d voxels)\n', ...
                volume.V_focal_mm3, volume.N_voxels);
        fprintf('Focal volume (isosurface):  %.4f mm^3\n', volume.V_isosurf_mm3);
        if isfinite(volume.V_isosurf_mm3)
            fprintf('Ratio isosurface/voxel:     %.4f\n', ...
                    volume.V_isosurf_mm3 / volume.V_focal_mm3);
        end
    end

    % ---------------- Optional visualisation ----------------
    if do_plot && ~isempty(fv_focal.vertices)
        % Padded index p maps to original index p-1, i.e. to physical
        % coordinate (p-2)*spacing.
        v = fv_focal.vertices;
        fv_plot = fv_focal;
        fv_plot.vertices = [ (v(:,1) - 2) * dy, ...
                             (v(:,2) - 2) * dx, ...
                             (v(:,3) - 2) * dz ];

        figure;
        patch(fv_plot, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        view(3); axis equal; camlight; lighting gouraud;
        xlabel('y (m)'); ylabel('x (m)'); zlabel('z (m)');
        title(sprintf('%.0f dB focal volume', opt.LevelDB));
        grid on;
    end
end

% --- Helper function ---
function V = isosurfaceVolume(vertices, faces)
    v1 = vertices(faces(:,1), :);
    v2 = vertices(faces(:,2), :);
    v3 = vertices(faces(:,3), :);
    V  = abs(sum(sum(v1 .* cross(v2, v3, 2), 2)) / 6);
end