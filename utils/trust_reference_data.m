function ref = trust_reference_data()
%TRUST_REFERENCE_DATA  Reference measurements and simulated focal gain for the
% TRUST array amplitude calibration.
%
% All pressures are returned in Pa. The source tables report MPa.
%
% Sources:
%   Tables 2.2 / 2.3 / 2.4 of the thesis, from TUV SUD report no. 713377467.
%   Simulated focal gain from the free-field acousticFieldPropagator sweep.
%
% The calibration metric is the effective pressure
%       p_eff = ( p_max + |p_min| ) / 2
% which is the symmetric part of the focal waveform and is therefore
% insensitive to the first-order nonlinear distortion visible in p_max.
%
%   author: Maylis Mueller

% -------------------------------------------------------------------------
% DEPTH SWEEP AT V_ref  (Table 2.2)
% -------------------------------------------------------------------------

ref.V_ref = 30;                                             % [V]

ref.z_hyd    = [ 20    30    40    50    60    70    80    90   100   120 ];   % [mm]
pmax_hyd_MPa = [1.34  1.40  1.39  1.24  1.13  1.08  1.01  0.94  0.94  0.83];
pmin_hyd_MPa = [1.29  1.33  1.25  1.14  1.02  0.93  0.86  0.80  0.77  0.69];

ref.peff_hyd = 0.5 * (pmax_hyd_MPa + abs(pmin_hyd_MPa)) * 1e6;   % [Pa]


% -------------------------------------------------------------------------
% VOLTAGE SWEEPS  (Tables 2.3 and 2.4)
% -------------------------------------------------------------------------
% Retained for verification of the proportionality assumption only. They are
% NOT used to build s(V); see trust_check_svlinearity.

ref.V_30    = [  10    20    30    40    50    60    70    80    90 ];        % [V]
pmax_30_MPa = [0.47  0.92  1.50  1.88  2.49  3.00  3.44  3.86  4.42];
pmin_30_MPa = [0.46  0.90  1.32  1.72  2.23  2.58  2.94  3.36  3.78];
ref.peff_30 = 0.5 * (pmax_30_MPa + abs(pmin_30_MPa)) * 1e6;      % [Pa]

ref.V_70    = [  10    20    30    40    50    60    70    80 ];              % [V]
pmax_70_MPa = [0.39  0.70  1.05  1.40  1.83  2.25  2.71  3.19];
pmin_70_MPa = [0.36  0.64  0.92  1.17  1.47  1.74  2.01  2.23];
ref.peff_70 = 0.5 * (pmax_70_MPa + abs(pmin_70_MPa)) * 1e6;      % [Pa]


% -------------------------------------------------------------------------
% SIMULATED FOCAL GAIN
% -------------------------------------------------------------------------
% Focal pressure produced by the calibrated source model per unit source
% amplitude, one forward run per commanded depth. In a linear medium without
% medium.BonA the focal waveform is symmetric, so p_eff_sim == p_max_sim.
% If p_min is recorded, replace peff_sim with 0.5*(p_max + |p_min|) as a check.

ref.p0_ref   = 1;                                                            % [Pa]
ref.z_sim    = [ 20    30    40    50    60    70    80    90   100   110   120 ];  % [mm]
ref.peff_sim = [13.90 17.12 17.62 16.93 15.85 14.72 13.67 12.73 11.91 11.19 10.57]; % [Pa]

ref.G = ref.peff_sim / ref.p0_ref;   % dimensionless focal gain of the model


% -------------------------------------------------------------------------
% VALIDITY RANGES
% -------------------------------------------------------------------------
% Voltage: p_eff/V is constant to within +/-2.5% over this range at both
% measured depths. Below 20 V the V=10 points are outliers in both sweeps;
% above 70 V the 70 mm response departs from proportionality.
%
% Depth: bounded by the measured depth sweep. No extrapolation.

ref.V_valid = [20  70];    % [V]
ref.z_valid = [20 120];    % [mm]

end
