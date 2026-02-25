function [data_p_rho] = get_p_rho_comment(LAI, SZA, SAA, va, angle, lidf, hc0, Height_c, Crowndeepth, gap_tot, gap_betw, gap_within, gap_S, gap_S_tot, gap_S_within, go_par, c, HotSpotPar, id, D)
% GET_P_RHO_COMMENT Calculates structural probabilities and 4-component GO gap fractions.
% 
% DESCRIPTION:
% This core function serves as the central hub for combining Geometric-Optical (GO) 
% theory and spectral invariants (p-theory). It calculates the four scene components 
% (sunlit/shaded crown, sunlit/shaded background) and integrates hemispherical 
% scattering probabilities for each viewing angle.
%
% INPUTS:
%   LAI         - Leaf Area Index of a single crown [m2/m2]
%   SZA, SAA    - Sun Zenith and Azimuth Angles [degrees]
%   va          - View angle matrix (N x 4)
%   angle       - Number of viewing angles
%   lidf        - Leaf Inclination Distribution Function vector
%   hc0         - Structural parameter / canopy height indicator
%   Height_c    - Crown center height [m]
%   Crowndeepth - Crown depth [m]
%   gap_* - Various gap fraction matrices (total, between-crown, within-crown)
%   go_par      - GO theory hotspot parameter
%   c           - Canopy-level clumping index
%   HotSpotPar  - Leaf-scale hotspot parameter
%   id          - Hemispherical interceptance (scene level)
%   D           - Ratio of diffuse to incoming irradiance
%
% OUTPUTS:
%   data_p_rho  - A 27 x N matrix containing integrated scattering probabilities, 
%                 transmittances, and 4-component gap fractions for each view angle.

% -------------------------------------------------------------------------
% Initialize output matrix (27 parameters for each viewing angle)
data_p_rho = zeros(27, angle);

for t = 1:angle
    
    tts = SZA;
    tto = va(t, 1);    % View zenith angle 
    psi = va(t, 2);    % View azimuth angle
    
    % Ensure relative azimuth angle is wrapped correctly and symmetric [0, 180]
    if psi > 180
        psi = psi - 360; 
    end
    psi = abs(psi);    
    psi = abs(psi - 360 * round(psi / 360));  
    
    % Set empirical CI parameters to 1 (assuming GO model handles macroscopic clumping)
    CIy1 = 1;
    CIy2 = 1;
    [CIs] = CIxy(CIy1, CIy2, tts);
    [CIo] = CIxy(CIy1, CIy2, tto);
    
    % ---------------------------------------------------------------------
    % 1. Directional gap fractions extraction
    gap_V_tot  = gap_tot(t, 3);
    gap_V_betw = gap_betw(t, 3);
    Ps_dir_go  = gap_S;         % Between-crown gap fraction in solar direction
    Pv_dir_go  = gap_V_betw;    % Between-crown gap fraction in view direction
    
    % ---------------------------------------------------------------------
    % 2. Calculate four GO components (Kc, Kt, Kg, Kz)
    
    % Sunlit background (Kg) with hotspot effect
    Kg = Ps_dir_go * Pv_dir_go + get_HSF_go(go_par, tts, SAA, tto, psi, Ps_dir_go, Pv_dir_go, Height_c);
    
    % Shaded background (Kz)
    Kz = Pv_dir_go - Kg;   % Note: Kg + Kz = F (total background seen)
    
    % Total crown seen
    Kct = 1 - Pv_dir_go;
    
    % Phase angle calculation for canopy shadowing
    delta_angle = cosd(tts) * cosd(tto) + sind(tts) * sind(tto) * cosd(psi - SAA); 
    phi_angle   = acosd(delta_angle);
    delta_val   = cosd(phi_angle .* (1 - sin(pi .* c ./ 2)));   
    
    % Sunlit crown (Kc)
    if ((hc0 - Crowndeepth) < Crowndeepth) && (tto > tts) && (SAA == psi)   
        % Continuous canopy scenario (no shaded crown visible)
        Kc = Kct; 
    else
        % Discrete canopy scenario
        Kc = 0.5 * (1 + delta_val) * Kct;  
    end
    
    % Shaded crown (Kt)
    Kt = Kct - Kc;
    
    % ---------------------------------------------------------------------
    % 3. Crown-level scattering and transmittance parameters
    Ps_dir_inKz = gap_S_within;     % Within-crown gap fraction in solar direction
    
    [Gs, Go, k, K, sob, sof] = PHASE(tts, tto, psi, lidf);
    
    % Sunlit/Shaded probabilities and scattering
    [kc, kg]       = sunshade_H(tts, tto, psi, Gs, Go, CIs, CIo, LAI, HotSpotPar);   
    [kc_kt, kg_kt] = sunshade_Kt_He(tts, tto, psi, Gs, Go, CIs, CIo, LAI);   
    
    % ---------------------------------------------------------------------
    % 4. Hemispherical integration and p-theory parameters
    
    i0 = 1 - gap_S_tot;
    i0 = D * id + (1 - D) * i0;
    iv = 1 - gap_V_tot;
    tv = 1 - iv;
    
    p         = 1 - id / LAI;
    rho2      = iv / 2 / LAI;
    rho_hemi2 = id / 2 / LAI;
    
    % Directional-hemispherical scattering (fixed solar angle)
    [sob_vsla, sof_vsla, kgd] = A_BRFv2_single_hemi(tts, CIs, CIy1, CIy2, LAI, lidf, hc0);  
    
    % Hemispherical-directional scattering (fixed view angle)
    [sob_vsla_dif, sof_vsla_dif, kg_dif] = A_BRFv2_single_dif(tto, CIo, CIy1, CIy2, LAI, lidf, hc0);  
    
    % Bi-hemispherical scattering
    [sob_vsla_hemi_dif, sof_vsla_hemi_dif, kgd_dif] = A_BRFv2_single_hemi_dif(CIy1, CIy2, LAI, lidf, hc0);  
    
    % ---------------------------------------------------------------------
    % 5. Pack variables into output vector
    data_p_rho(:, t) = [i0, id, p, rho_hemi2, sob_vsla, sof_vsla, rho2, tv, ...
                        sob, sof, kc, kc_kt, kg, kg_kt, K, kgd, ...
                        sob_vsla_dif, sof_vsla_dif, kg_dif, ...
                        sob_vsla_hemi_dif, sof_vsla_hemi_dif, kgd_dif, ...
                        Ps_dir_inKz, Kg, Kc, Kt, Kz]';
    
end
end