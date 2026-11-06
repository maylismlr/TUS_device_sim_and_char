function [karray, positions] = build_karray(row_counts, element_height, element_width, element_pitch, z_pos)
%BUILD_KARRAY  Geometry of the 256-element TRUST array as a kWaveArray object.
%
%   [KARRAY, POSITIONS] = BUILD_KARRAY(ROW_COUNTS, ELEMENT_SIZE, ELEMENT_PITCH,
%   Z_POS) places square elements on a regular grid, row by row, centred on the
%   origin in x and y, with the array face at Z_POS.
%
%   This function exists so that there is exactly ONE definition of the array
%   geometry in the project. create_source and fit_element_width both call it.
%   Previously each built its own kWaveArray, which meant the element size
%   fitted by one could silently stop matching the array simulated by the other.
%
%   No kgrid argument: a kWaveArray stores continuous coordinates, not grid
%   indices -- that is the whole point of the off-grid source method. The grid
%   only enters later, when the caller asks for grid weights or a binary mask.
%   This also means the positions can be queried before a grid exists, which is
%   how fit_element_width sizes its domain.
%
%   Inputs:
%     row_counts     number of elements in each row, e.g.
%                    [6 10 12 14 16 16 18 18 18 18 18 18 16 16 14 12 10 6]
%     element_size   side length of the square elements [m]
%     element_pitch  centre-to-centre spacing [m]
%     z_pos          axial position of the array face [m]
%
%   Outputs:
%     karray         kWaveArray with sum(row_counts) rectangular elements
%     positions      [N x 3] element centre coordinates [m]
%
%   Note on element_size: this is a free parameter, not a physical dimension.
%   The physical elements are 2.58 mm squares on a 3.08 mm pitch, but the
%   EFFECTIVE radiating area is larger, because the element vibration couples
%   into the surrounding structure and the polyurethane coupling layer smears
%   the source across the kerf. fit_element_width determines the effective value
%   from measurements.
%
%   author: Maylis Mueller

% Rows are centred on the origin. The (num_rows + 1)/2 form handles both
% parities: with an odd number of rows the centre falls on the middle row, with
% an even number it falls halfway between the two middle rows.
num_rows = numel(row_counts);
y_center = (num_rows + 1)/2;

% 'BLITolerance' and 'UpsamplingRate' control the off-grid source representation
% (Wise et al. 2019): each element is expanded as a band-limited interpolant
% rather than staircased onto grid points, which is what makes a 2.58 mm element
% meaningful on a 0.5 mm grid at all. A lower BLITolerance gives a wider support
% and a more accurate element, at the cost of speed and memory.
karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

positions = zeros(sum(row_counts), 3);
el_idx = 1;

for i = 1:num_rows

    n_el = row_counts(i);

    % y position of this row
    y_pos = (i - y_center) * element_pitch;

    % elements within the row are centred in x by the same rule
    x_center = (n_el + 1)/2;

    for j = 1:n_el

        x_pos = (j - x_center) * element_pitch;

        positions(el_idx, :) = [x_pos, y_pos, z_pos];

        % addRectElement(centre, Lx, Ly, rotation). The [0 0 0] rotation leaves
        % the element in the xy plane, i.e. radiating along +z.
        karray.addRectElement([x_pos, y_pos, z_pos], ...
                              element_height, element_width, [0 0 0]);

        el_idx = el_idx + 1;

    end
end

end