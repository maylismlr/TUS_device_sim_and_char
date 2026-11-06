function focal_mask = get_focal_mask(p_field, focus_idx, p_fraction)
%GET_FOCAL_MASK
% Returns only the connected thresholded component containing the
% supplied focal peak.
%
% For -6 dB pressure:
%       p_fraction = 0.5
% For -3 dB pressure:
%       p_fraction = 0.707
%
%   author: Maylis Mueller
    ix = focus_idx(1);
    iy = focus_idx(2);
    iz = focus_idx(3);

    p_peak = p_field(ix, iy, iz);

    threshold = p_fraction * p_peak;

    mask = p_field >= threshold;

    CC = bwconncomp(mask, 26);

    focus_linear = sub2ind( ...
        size(p_field), ...
        ix, iy, iz);

    focal_mask = false(size(mask));

    for c = 1:CC.NumObjects

        if any(CC.PixelIdxList{c} == focus_linear)

            focal_mask(CC.PixelIdxList{c}) = true;
            return;

        end

    end

    error('Could not identify connected focal component.');

end

