function [amp_in, phase_in] = build_afp_source_bli( ...
    positions, element_size, source_amp, ...
    x_vec, y_vec, z_vec, focal_depth, f0, c0)
%BUILD_AFP_SOURCE_BLI
% Source amplitude and phase matrices for acousticFieldPropagator, with the
% elements represented as band-limited interpolants instead of being
% rasterised onto the grid.
%
%
% WHY
% The rasterising version tests each grid point for membership:
%
%   x_ind = find(x_vec >= x_element - h & x_vec <= x_element + h);
%
% That is a binary in/out decision, so the represented element width can only
% change in whole grid steps. At dx = 0.25 mm a 3.0 mm and a 3.2 mm element
% select the same 13 cells whenever the element centre lands on a grid point,
% and the two candidates produce a byte-identical source. Sweeping element
% size in 0.1 mm steps then asks a question the representation cannot answer.
%
% kWaveArray instead computes the element's spatial Fourier content and lays
% down fractional (and negative) weights that reproduce that band-limited
% shape. Two sizes differing by less than one grid step give genuinely
% different weights, so element size is resolved below dx without refining
% the grid (Wise et al., JASA 146(1), 2019).
%
% A second, quieter benefit: element centres sit at multiples of the pitch,
% which is not generally a multiple of dx, so under rasterisation different
% elements round to different widths. That area jitter across the aperture is
% a numerical artefact. The band-limited representation does not have it.
%
% WHY 2D
% The source occupies a single z plane, and a planar element's 3D band-limited
% representation factorises into a 2D interpolant in xy times a delta in z.
% So the array is built on a 2D grid and the result written into plane 1.
% Building it on the full 3D grid would compute the same thing 581 times over.
%
% NOTE ON SCALE
% kWaveArray weights carry their own normalisation, so the absolute source
% strength differs from the rasterised version by a constant. That constant
% cancels in the normalised depth curve. It does NOT cancel in the absolute
% amplitude calibration, which must be redone against whichever source
% definition is finally used.
%
% author: Maylis Mueller

Nx = numel(x_vec);
Ny = numel(y_vec);
Nz = numel(z_vec);
dx = x_vec(2) - x_vec(1);

n_elements = size(positions, 1);

% -------------------------------------------------------------------------
% Element weights.
%
% These depend only on the geometry, not on the focal depth, but the caller
% invokes this function once per focus. Caching them means the band-limited
% expansion is computed once per candidate rather than once per candidate per
% depth -- a tenfold saving on the expensive part.
% -------------------------------------------------------------------------
persistent cached_weights cached_positions cached_size cached_grid

grid_signature = [Nx, Ny, dx];

rebuild = isempty(cached_weights) ...
       || ~isequal(cached_positions, positions) ...
       || ~isequal(cached_size, element_size) ...
       || ~isequal(cached_grid, grid_signature);

if rebuild

    kgrid = kWaveGrid(Nx, dx, Ny, dx);

    % Sanity check: the caller's axes must match the kWaveGrid axes, or the
    % element positions will be written to the wrong grid points. Both use
    % ((0:N-1) - floor(N/2))*dx for odd N, so they should agree exactly.
    if max(abs(kgrid.x_vec(:) - x_vec(:))) > 1e-12 || ...
       max(abs(kgrid.y_vec(:) - y_vec(:))) > 1e-12
        error('build_afp_source_bli:axisMismatch', ...
              ['kWaveGrid axes do not match x_vec/y_vec. Use an odd number ' ...
               'of transverse grid points so that zero falls on a grid point.']);
    end

    karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

    for element_ind = 1:n_elements
        % 2D addRectElement takes [x y], side lengths, and a scalar rotation.
        karray.addRectElement( ...
            [positions(element_ind, 1), positions(element_ind, 2)], ...
            element_size, element_size, 0);
    end

    % Store sparsely. Each element's support is local, so holding 256 full
    % Nx-by-Ny matrices would be wasteful.
    cached_weights = cell(n_elements, 1);

    for element_ind = 1:n_elements
        w  = karray.getElementGridWeights(kgrid, element_ind);
        ix = find(w ~= 0);
        cached_weights{element_ind} = struct('idx', ix, 'val', w(ix));
    end

    cached_positions = positions;
    cached_size      = element_size;
    cached_grid      = grid_signature;

    fprintf('    [BLI] rebuilt weights for %.3f mm elements (%d nonzero points/element)\n', ...
            element_size*1e3, round(mean(cellfun(@(c) numel(c.idx), cached_weights))));
end

% -------------------------------------------------------------------------
% Focusing phases
% -------------------------------------------------------------------------
k0 = 2*pi*f0/c0;

distance_to_focus = sqrt( ...
    positions(:, 1).^2 + ...
    positions(:, 2).^2 + ...
    focal_depth.^2);

% Only relative phase matters; subtracting the minimum keeps the numbers small.
element_phase = +k0 * (distance_to_focus - min(distance_to_focus));

% -------------------------------------------------------------------------
% Superpose. Weights can be negative (the band-limited interpolant rings), so
% accumulate complex and take abs/angle at the end -- a negative weight is
% simply |w| with a phase offset of pi, which angle() handles correctly.
% Overlapping elements add, as before.
% -------------------------------------------------------------------------
source_complex_2D = complex(zeros(Nx, Ny));

for element_ind = 1:n_elements
    w = cached_weights{element_ind};
    source_complex_2D(w.idx) = source_complex_2D(w.idx) ...
        + source_amp * w.val * exp(1i * element_phase(element_ind));
end

amp_in   = zeros(Nx, Ny, Nz);
phase_in = zeros(Nx, Ny, Nz);

amp_in(:, :, 1)   = abs(source_complex_2D);
phase_in(:, :, 1) = angle(source_complex_2D);

end