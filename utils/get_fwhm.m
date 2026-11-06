function [fwhm, z_left, z_right] = get_fwhm(z, p, peak_ind)
%   author: Maylis Mueller

    z = z(:);
    p = p(:);

    half_max = 0.5 * p(peak_ind);

    % Left crossing
    i = peak_ind;

    while i > 1 && p(i) >= half_max
        i = i - 1;
    end

    if i == 1 && p(i) >= half_max
        z_left = NaN;
    else
        z_left = interp1( ...
            p(i:i+1), ...
            z(i:i+1), ...
            half_max, ...
            'linear');
    end

    % Right crossing
    i = peak_ind;

    while i < numel(p) && p(i) >= half_max
        i = i + 1;
    end

    if i == numel(p) && p(i) >= half_max
        z_right = NaN;
    else
        z_right = interp1( ...
            p(i-1:i), ...
            z(i-1:i), ...
            half_max, ...
            'linear');
    end

    fwhm = z_right - z_left;

end