function [p0, info] = trust_p0(V_input, z_focus, ref)
%TRUST_P0  Source amplitude for a k-Wave TRUST simulation.
%
%   p0 = TRUST_P0(V, z) returns the source amplitude [Pa] required so that the
%   simulated effective focal pressure matches that measured on the TRUST
%   device driven at V volts with a commanded focal depth of z mm.
%
%   [p0, info] = TRUST_P0(...) also returns the two factors separately.
%
%   The calibration factorises as
%
%       p0(V,z) = s(V) * k(z)
%
%   s(V) = V / V_ref                                          [dimensionless]
%       Voltage scale, a property of the drive electronics. It carries no
%       depth argument: the beamformer sets element delays and not element
%       amplitudes, so the volts-to-surface-pressure map is taken to be
%       independent of where the array is steered. The response is
%       proportional through the origin, so the slope cancels in the ratio and
%       no coefficient is fitted.
%
%   k(z) = p_eff_meas(z, V_ref) / G_sim(z)                              [Pa]
%       Source amplitude required at the reference drive, where G_sim is the
%       dimensionless model focal gain. Since
%       k(z) = f(V_ref) * G_dev(z) / G_sim(z), it would be constant with depth
%       if the model reproduced the device's focal-gain variation; its
%       variation therefore quantifies the residual model error. Indexed by the
%       NOMINAL COMMANDED depth, not by any beamformer-corrected depth, so that
%       it composes with the geometric correction without double counting.
%
%   The arbitrary reference amplitude p0_ref used to obtain G_sim cancels and
%   does not appear here.
%
%   INPUTS
%       V_input  Excitation voltage [V], scalar or array
%       z_focus  Commanded focal depth [mm], scalar or array (same size, or
%                one of the two may be scalar)
%       ref      (optional) struct from trust_reference_data
%
%   OUTPUT
%       p0       Source amplitude [Pa]
%       info     struct with fields s, k, k_nodes, G, p_target
%
%   author: Maylis Mueller

if nargin < 3 || isempty(ref)
    ref = trust_reference_data();
end

% Allow scalar expansion of either argument
if isscalar(V_input) && ~isscalar(z_focus)
    V_input = repmat(V_input, size(z_focus));
elseif isscalar(z_focus) && ~isscalar(V_input)
    z_focus = repmat(z_focus, size(V_input));
end
if ~isequal(size(V_input), size(z_focus))
    error('trust_p0:sizeMismatch', ...
          'V_input and z_focus must be the same size, or one must be scalar.');
end


% -------------------------------------------------------------------------
% RANGE CHECKS  (hard error on depth, warning on voltage)
% -------------------------------------------------------------------------

if any(z_focus(:) < ref.z_valid(1)) || any(z_focus(:) > ref.z_valid(2))
    error('trust_p0:depthOutOfRange', ...
          ['Commanded depth outside the measured range [%g %g] mm. ' ...
           'k(z) is not defined there and must not be extrapolated.'], ...
          ref.z_valid(1), ref.z_valid(2));
end

if any(V_input(:) < ref.V_valid(1)) || any(V_input(:) > ref.V_valid(2))
    warning('trust_p0:voltageOutOfRange', ...
            ['Voltage outside the range [%g %g] V over which p_eff is ' ...
             'proportional to V. s(V) = V/V_ref may not hold.'], ...
            ref.V_valid(1), ref.V_valid(2));
end


% -------------------------------------------------------------------------
% 1. VOLTAGE SCALE  s(V)
% -------------------------------------------------------------------------

s = V_input / ref.V_ref;


% -------------------------------------------------------------------------
% 2. DEPTH SCALE  k(z)
% -------------------------------------------------------------------------
% k is formed on the measured depth grid, then interpolated. Interpolating k
% itself rather than measurement and gain separately avoids two independent
% interpolation errors entering the quotient. Note that every depth in z_hyd
% is also present in z_sim, so the first interp1 returns the nodes exactly and
% introduces no error of its own.

G_at_hyd = interp1(ref.z_sim, ref.G, ref.z_hyd, 'pchip');   % model gain on the
                                                            % measured grid
k_nodes = ref.peff_hyd ./ G_at_hyd;   % [Pa]

k = interp1(ref.z_hyd, k_nodes, z_focus, 'pchip');



% -------------------------------------------------------------------------
% 3. SOURCE AMPLITUDE
% -------------------------------------------------------------------------

p0 = s .* k;                          % [Pa]

if nargout > 1
    info.s        = s;
    info.k        = k;
    info.k_nodes  = k_nodes;
    info.G        = interp1(ref.z_sim, ref.G, z_focus, 'pchip');
    info.p_target = s .* interp1(ref.z_hyd, ref.peff_hyd, z_focus, 'pchip');
end

end
