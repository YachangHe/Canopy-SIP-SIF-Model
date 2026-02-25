function [sob_vsla, sof_vsla, kgd_dif] = A_BRFv2_single_hemi_dif(CIy1, CIy2, LAI, lidf, hc0)
% A_BRFV2_SINGLE_HEMI_DIF Calculates bi-hemispherical single scattering and gap probabilities.
% 
% DESCRIPTION:
% This function performs a 4D numerical integration over both the incident (solar) 
% hemisphere and the scattered (viewing) hemisphere using an 8-point Gaussian 
% quadrature method. It computes the bi-hemispherical (diffuse-diffuse) integration 
% of scattering and transmittance.
%
% INPUTS:
%   CIy1 - Parameter 1 for the angular dependence of the clumping index
%   CIy2 - Parameter 2 for the angular dependence of the clumping index
%   LAI  - Leaf Area Index [m2/m2]
%   lidf - Leaf Inclination Distribution Function vector
%   hc0  - Hotspot parameter or structural parameter related to canopy height
%
% OUTPUTS:
%   sob_vsla - Integrated backward single scattering contribution (bi-hemispherical)
%   sof_vsla - Integrated forward single scattering contribution (bi-hemispherical)
%   kgd_dif  - Integrated bi-hemispherical gap probability (diffuse-diffuse transmittance)

% -------------------------------------------------------------------------
% 8-point Gaussian quadrature nodes (xx) and weights (ww)
xx = [0.9602898565 -0.9602898565 0.7966664774 -0.7966664774 0.5255324099 -0.5255324099 0.1834346425 -0.1834346425];
ww = [0.1012285363 0.1012285363 0.2223810345 0.2223810345 0.3137066459 0.3137066459 0.3626837834 0.3626837834];

% 1. Integration limits for INCIDENT Zenith Angle (theta_in -> mL)
upperlimit_mL = pi/2.0;
lowerlimit_mL = 0.0;
conv1_mL = (upperlimit_mL - lowerlimit_mL) / 2.0;
conv2_mL = (upperlimit_mL + lowerlimit_mL) / 2.0;

% 2. Integration limits for INCIDENT Azimuth Angle (phi_in -> nL)
upperlimit_nL = 2.0 * pi;
lowerlimit_nL = 0.0;
conv1_nL = (upperlimit_nL - lowerlimit_nL) / 2.0;
conv2_nL = (upperlimit_nL + lowerlimit_nL) / 2.0;

% 3. Integration limits for SCATTERED Zenith Angle (theta_out -> tL)
upperlimit_tL = pi/2.0;
lowerlimit_tL = 0.0;
conv1_tL = (upperlimit_tL - lowerlimit_tL) / 2.0;
conv2_tL = (upperlimit_tL + lowerlimit_tL) / 2.0;

% 4. Integration limits for SCATTERED Azimuth Angle (phi_out -> pL)
upperlimit_pL = 2.0 * pi;
lowerlimit_pL = 0.0;
conv1_pL = (upperlimit_pL - lowerlimit_pL) / 2.0;
conv2_pL = (upperlimit_pL + lowerlimit_pL) / 2.0;

% Initialize outermost accumulators
sum_mL   = 0.0;
sum_mL_f = 0.0;
sum_mL_g = 0.0;

% --- Loop 1: Incident Zenith Angle (m) ---
for m = 1:8
    
    neword_mL = conv1_mL * xx(m) + conv2_mL;
    mu_mL     = cos(neword_mL);
    sin_mL    = sin(neword_mL);
    
    sum_nL   = 0.0;
    sum_nL_f = 0.0;
    sum_nL_g = 0.0;
    
    % --- Loop 2: Incident Azimuth Angle (n) ---
    for n = 1:8
        neword_nL = conv1_nL * xx(n) + conv2_nL;
        
        sum_tL   = 0.0;
        sum_tL_f = 0.0;
        sum_tL_g = 0.0;
        
        % --- Loop 3: Scattered Zenith Angle (i) ---
        for i = 1:8
            
            neword_tL = conv1_tL * xx(i) + conv2_tL;
            mu_tL     = cos(neword_tL);
            sin_tL    = sin(neword_tL);
            
            sum_pL   = 0.0;
            sum_pL_f = 0.0;
            sum_pL_g = 0.0;
            
            % --- Loop 4: Scattered Azimuth Angle (j) ---
            for j = 1:8
                
                neword_pL = conv1_pL * xx(j) + conv2_pL;
                
                % Convert to degrees for physical subroutines
                tts = neword_mL * 180 / pi;
                tto = neword_tL * 180 / pi;
                
                % Calculate relative azimuth angle
                psi = abs(neword_nL * 180 / pi - neword_pL * 180 / pi);
                psi = abs(psi - 360 * round(psi / 360));
                
                % Calculate phase functions and structural parameters
                [Gs, Go, k, K, sob, sof] = PHASE(tts, tto, psi, lidf);
                
                [CIs] = CIxy(CIy1, CIy2, tts);
                [CIo] = CIxy(CIy1, CIy2, tto);
                
                [kca, kga] = sunshade(tts, tto, psi, Gs, Go, CIs, CIo, LAI, hc0);
                
                % Accumulate innermost integral
                sum_pL   = sum_pL   + ww(j) * sob .* kca / K / pi;
                sum_pL_f = sum_pL_f + ww(j) * sof .* kca / K / pi;
                sum_pL_g = sum_pL_g + ww(j) * kga / pi;
                
            end
            
            % Apply Loop 4 conversions and accumulate into Loop 3
            sum_pL   = sum_pL * conv1_pL;
            sum_tL   = sum_tL + ww(i) * mu_tL * sin_tL * sum_pL;
            
            sum_pL_f = sum_pL_f * conv1_pL;
            sum_tL_f = sum_tL_f + ww(i) * mu_tL * sin_tL * sum_pL_f;
            
            sum_pL_g = sum_pL_g * conv1_pL;
            sum_tL_g = sum_tL_g + ww(i) * mu_tL * sin_tL * sum_pL_g;
        end
        
        % Apply Loop 3 conversions and accumulate into Loop 2
        sum_tL   = sum_tL * conv1_tL;
        sum_nL   = sum_nL + ww(n) * sum_tL / pi;
        
        sum_tL_f = sum_tL_f * conv1_tL;
        sum_nL_f = sum_nL_f + ww(n) * sum_tL_f / pi;
        
        sum_tL_g = sum_tL_g * conv1_tL;
        sum_nL_g = sum_nL_g + ww(n) * sum_tL_g / pi;
    end
    
    % Apply Loop 2 conversions and accumulate into Loop 1
    sum_nL   = sum_nL * conv1_nL;
    sum_mL   = sum_mL + ww(m) * mu_mL * sin_mL * sum_nL;
    
    sum_nL_f = sum_nL_f * conv1_nL;
    sum_mL_f = sum_mL_f + ww(m) * mu_mL * sin_mL * sum_nL_f;
    
    sum_nL_g = sum_nL_g * conv1_nL;
    sum_mL_g = sum_mL_g + ww(m) * mu_mL * sin_mL * sum_nL_g;
end

% Apply Loop 1 conversions to finalize integration
sum_mL   = sum_mL * conv1_mL;
sum_mL_f = sum_mL_f * conv1_mL;
sum_mL_g = sum_mL_g * conv1_mL;

% Assign to output variables
sob_vsla = sum_mL;
sof_vsla = sum_mL_f;
kgd_dif  = sum_mL_g;
end