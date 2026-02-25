function [bg] = get_HSF_go(par, SZA, SAA, VZA, VAA, Ps_dir_go, Pv_dir_go, z)
% GET_HSF_GO Calculates the bidirectional between-crown gap probability covariance.
% 
% DESCRIPTION:
% This routine calculates the hotspot enhancement term (covariance) for the 
% bidirectional between-crown gap probability based on Geometric-Optical (GO) theory.
% 
% The total bidirectional gap probability is defined as:
% KG = Ps * Pv + bg, where bg = Y * sqrt(Ps * Pv * (1 - Ps) * (1 - Pv))
% This function specifically returns the covariance term 'bg' (Equation 7 in the paper).
%
% REFERENCE:
% This function is adapted from the PATH_RT model:
% https://doi.org/10.1016/j.rse.2023.113985
%
% INPUTS:
%   par       - Crown size or hotspot parameter [m]
%   SZA       - Sun Zenith Angle [degrees]
%   SAA       - Sun Azimuth Angle [degrees]
%   VZA       - View Zenith Angle [degrees]
%   VAA       - View Azimuth Angle [degrees]
%   Ps_dir_go - Directional between-crown gap probability in solar direction
%   Pv_dir_go - Directional between-crown gap probability in view direction
%   z         - Crown center height [m] (Height - lmax + 1/2*Crowndeepth)
%
% OUTPUTS:
%   bg        - The covariance term of the bidirectional gap probability

% -------------------------------------------------------------------------
% Convert angles from degrees to radians
SZA_rad = deg2rad(SZA);
SAA_rad = deg2rad(SAA);
VZA_rad = deg2rad(VZA);
VAA_rad = deg2rad(VAA);

% Calculate the covariance factor (f1)
f1 = sqrt(Ps_dir_go * Pv_dir_go * (1 - Ps_dir_go) * (1 - Pv_dir_go));

% Calculate the cosine of the scattering angle (phase angle)
cosgamma = cos(SZA_rad) * cos(VZA_rad) + sin(SZA_rad) * sin(VZA_rad) * cos(VAA_rad - SAA_rad);

% Calculate the angular distance measure (delta) based on Equation (9)
delta = sqrt(1 / (cos(SZA_rad)^2) + 1 / (cos(VZA_rad)^2) - 2 * cosgamma / (cos(SZA_rad) * cos(VZA_rad)));

% Apply a small threshold to delta to ensure numerical stability (prevent exactly zero)
if delta < 0.00001
    delta = 0.00001;
end

% Calculate the overlap function (Y)
Y = exp(-delta / par * z);

% Calculate the final hotspot enhancement term
bg = f1 * Y;

end