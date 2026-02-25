function [sob_vsla, sof_vsla, kgd] = A_BRFv2_single_hemi(tts, CIs, CIy1, CIy2, LAI, lidf, hc0)
% A_BRFV2_SINGLE_HEMI Calculates directional-hemispherical single scattering and gap probabilities.
% 
% DESCRIPTION:
% This function performs a double numerical integration over the upper hemisphere 
% (zenith angle: 0 to pi/2, azimuth angle: 0 to 2*pi) using an 8-point Gaussian 
% quadrature method. It computes the hemispherical integration of scattering and 
% transmittance for a specific solar incidence angle.
%
% INPUTS:
%   tts  - Solar zenith angle [degrees]
%   CIs  - Clumping index at the solar direction
%   CIy1 - Parameter 1 for the angular dependence of the clumping index
%   CIy2 - Parameter 2 for the angular dependence of the clumping index
%   LAI  - Leaf Area Index [m2/m2]
%   lidf - Leaf Inclination Distribution Function vector
%   hc0  - Hotspot parameter or structural parameter related to canopy height
%
% OUTPUTS:
%   sob_vsla - Integrated backward single scattering contribution
%   sof_vsla - Integrated forward single scattering contribution
%   kgd      - Integrated gap probability (directional-hemispherical transmittance)

% -------------------------------------------------------------------------
% 8-point Gaussian quadrature nodes (xx) and weights (ww)
xx = [0.9602898565 -0.9602898565 0.7966664774 -0.7966664774 0.5255324099 -0.5255324099 0.1834346425 -0.1834346425];
ww = [0.1012285363 0.1012285363 0.2223810345 0.2223810345 0.3137066459 0.3137066459 0.3626837834 0.3626837834];

% Define limits of integration and the conversion factors for integration
% over thetaL (Scattered Zenith Angle integration limits: 0 to pi/2)
upperlimit_tL = pi/2.0;
lowerlimit_tL = 0.0;
conv1_tL = (upperlimit_tL - lowerlimit_tL) / 2.0;
conv2_tL = (upperlimit_tL + lowerlimit_tL) / 2.0;

% Define limits of integration and the conversion factors for integration
% over phiL (Scattered Azimuth Angle integration limits: 0 to 2*pi)
upperlimit_pL = 2.0 * pi;
lowerlimit_pL = 0.0;
conv1_pL = (upperlimit_pL - lowerlimit_pL) / 2.0;
conv2_pL = (upperlimit_pL + lowerlimit_pL) / 2.0;

% Initialize integration accumulators
sum_tL   = 0;
sum_tL_f = 0;
sum_tL_g = 0;

% Outer loop: Integration over zenith angle (thetaL)
for i = 1:8
    
    neword_tL = conv1_tL * xx(i) + conv2_tL;
    mu_tL     = cos(neword_tL);
    sin_tL    = sin(neword_tL);
    
    sum_pL   = 0;
    sum_pL_f = 0;
    sum_pL_g = 0;
    
    % Inner loop: Integration over azimuth angle (phiL)
    for j = 1:8
        
        neword_pL = conv1_pL * xx(j) + conv2_pL;
        
        tta  = neword_tL * 180 / pi;   % Convert scattered zenith angle to degrees
        psia = neword_pL * 180 / pi;   % Convert relative azimuth angle to degrees
        
        % Calculate scattering phase functions and geometrical parameters
        [Gs, Ga, k, Ka, sob, sof] = PHASE(tts, tta, psia, lidf);      
        
        % Calculate angular clumping index
        [CIa] = CIxy(CIy1, CIy2, tta);
        
        % Calculate sunlit and shaded probabilities
        [kca, kga] = sunshade(tts, tta, psia, Gs, Ga, CIs, CIa, LAI, hc0);     
        
        % Accumulate inner integral (weighted by Gaussian weights)
        sum_pL   = sum_pL   + ww(j) * sob .* kca / Ka / pi;
        sum_pL_f = sum_pL_f + ww(j) * sof .* kca / Ka / pi;
        sum_pL_g = sum_pL_g + ww(j) * kga / pi;
    end
    
    sum_pL = sum_pL * conv1_pL;
    sum_tL = sum_tL + ww(i) * mu_tL * sin_tL * sum_pL;
    
    sum_pL_f = sum_pL_f * conv1_pL;
    sum_tL_f = sum_tL_f + ww(i) * mu_tL * sin_tL * sum_pL_f;
    
    sum_pL_g = sum_pL_g * conv1_pL;
    sum_tL_g = sum_tL_g + ww(i) * mu_tL * sin_tL * sum_pL_g;
    
end

% Finalize integration by applying outer conversion factors
sum_tL   = sum_tL * conv1_tL;
sum_tL_f = sum_tL_f * conv1_tL;
sum_tL_g = sum_tL_g * conv1_tL;

% Assign to output variables
sob_vsla = sum_tL;
sof_vsla = sum_tL_f;
kgd      = sum_tL_g;
end