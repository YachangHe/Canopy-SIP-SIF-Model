function [Qfdir_all, Qfyld_all, Qapar_all] = get_rta_withsif_comment(nb, nf, ne, iwlfi, iwlfo, angle, Qins, Qind, rho_l, tau_l, rs, MfI, MbI, data_p_rho)
% GET_RTA_WITHSIF_COMMENT Calculates radiative transfer with SIF emission.
%
% DESCRIPTION:
% This function simulates the radiative transfer within the canopy including 
% the excitation and scattering of Solar-Induced chlorophyll Fluorescence (SIF). 
% It solves the radiative transfer equation by separating it into the Black 
% Soil (BS) problem and the Soil (S) problem. 
%
% INPUTS:
%   nb, nf, ne - Number of bands for total, fluorescence, and excitation spectra
%   iwlfi, iwlfo - Indices for excitation and observation bands
%   angle      - Number of view angles
%   Qins, Qind - Incident direct and diffuse radiation spectra
%   rho_l, tau_l - Leaf-level reflectance and transmittance spectra
%   rs         - Soil reflectance spectrum
%   MfI, MbI   - Fluorescence excitation matrices (forward and backward)
%   data_p_rho - Structural parameters from p-theory (27 x angle matrix)
%
% OUTPUTS:
%   Qfdir_all  - Total directional exitance including SIF (Top of Canopy)
%   Qfyld_all  - Total SIF yield within the canopy
%   Qapar_all  - Total Absorbed Photosynthetically Active Radiation (APAR)

% -------------------------------------------------------------------------
% Extract structural parameters from p-theory
i0 = data_p_rho(1, 1);       % Interceptance for solar direction
id = data_p_rho(2, 1);       % Hemispherical interceptance
p  = data_p_rho(3, 1);       % Recollision probability
rho2 = data_p_rho(4, 1);     % Hemispherical escape probability
sob_vsla = data_p_rho(5, 1); 
sof_vsla = data_p_rho(6, 1); 

rho = data_p_rho(7, :);      % Directional escape probability
tv  = data_p_rho(8, :);      % Transmittance in view direction
sob = data_p_rho(9, :);      % Phase function (backward)
sof = data_p_rho(10, :);     % Phase function (forward)
kc  = data_p_rho(11, :);     % Probability of viewing sunlit foliage
kc_kt = data_p_rho(12, :);
kg  = data_p_rho(13, :);
kg_kt = data_p_rho(14, :);
K   = data_p_rho(15, :);     % Extinction coefficient (observer)
kgd = data_p_rho(16, :); 

sob_vsla_dif = data_p_rho(17, :); 
sof_vsla_dif = data_p_rho(18, :); 
kg_dif       = data_p_rho(19, :); 
sob_vsla_hemi_dif = data_p_rho(20, 1); 
sof_vsla_hemi_dif = data_p_rho(21, 1);
kgd_dif      = data_p_rho(22, :); 
Ps_dir_inKz  = data_p_rho(23, :); 

Kg = data_p_rho(24, :); 
Kc = data_p_rho(25, :); 
Kt = data_p_rho(26, :);
Kz = data_p_rho(27, :); 

% Expand 4-component gap fractions across all spectral bands
Kg = repmat(Kg, nb, 1);
Kc = repmat(Kc, nb, 1);
Kt = repmat(Kt, nb, 1);
Kz = repmat(Kz, nb, 1);
Ps_dir_inKz = repmat(Ps_dir_inKz, nb, 1);

t0 = 1 - i0;
td = 1 - id; 
wleaf = rho_l + tau_l;  % Leaf single scattering albedo

% Total fluorescence excitation matrix
Mf = MfI + MbI;

% -------------------------------------------------------------------------
% 1. Black Soil (BS) Problem - Direct Illumination
% -------------------------------------------------------------------------
Qfdir = zeros(nb, angle, 11); 
Qfhemi = zeros(nb, 11);    
Qfdir_v1 = zeros(nb, angle, 11);
Qfdir_kt = zeros(nb, angle, 11);
Qfhemi_v1 = zeros(nb, 11);
Qfhemi_kt = zeros(nb, 11);

Qapar = zeros(nb, 11);   
Qdown = zeros(nb, 11);  
Qsig  = zeros(nb, 11);     
Qfyld = zeros(nf, 11);       

% First collision with incident direct light
Qsig(:, 1) = Qins * i0;  

for i = 1:11
    Qapar(:, i) = Qsig(:, i) .* (1 - wleaf);  
    MQ = Mf * Qsig(iwlfi, i);   % SIF excitation after collision
    
    if i == 1   % Single Scattering
        % Sunlit Crown
        Qfdir_v1(:, :, i) = Qins .* rho_l * (sob .* kc ./ K) + Qins .* tau_l * (sof .* kc ./ K);    
        Qfdir_v1(iwlfo, :, i) = Qfdir_v1(iwlfo, :, i) + MbI * Qins(iwlfi, 1) * (sob .* kc ./ K) + MfI * Qins(iwlfi, 1) * (sof .* kc ./ K); 
        Qfdir_v1(iwlfo, :, i) = Kc(iwlfo, :) .* Qfdir_v1(iwlfo, :, i); 
        
        % Shaded Crown
        Qfdir_kt(:, :, i) = Qins .* rho_l * (sob .* kc_kt ./ K) + Qins .* tau_l * (sof .* kc_kt ./ K);     
        Qfdir_kt(iwlfo, :, i) = Qfdir_kt(iwlfo, :, i) + MbI * Qins(iwlfi, 1) * (sob .* kc_kt ./ K) + MfI * Qins(iwlfi, 1) * (sof .* kc_kt ./ K); 
        Qfdir_kt(iwlfo, :, i) = Kt(iwlfo, :) .* sqrt(Ps_dir_inKz(iwlfo, :)) .* Qfdir_kt(iwlfo, :, i);  
        
        Qfdir(iwlfo, :, i) = Qfdir_v1(iwlfo, :, i) + Qfdir_kt(iwlfo, :, i);
        
        % Hemispherical components
        Qfhemi_v1(:, i) = Qins .* rho_l * sob_vsla + Qins .* tau_l * sof_vsla;
        Qfhemi_v1(iwlfo, i) = Qfhemi_v1(iwlfo, i) + MbI * Qins(iwlfi, 1) * sob_vsla + MfI * Qins(iwlfi, 1) * sof_vsla; 
        
        Qfhemi_kt(:, i) = Qins .* rho_l * sob_vsla + Qins .* tau_l * sof_vsla;
        Qfhemi_kt(iwlfo, i) = Qfhemi_kt(iwlfo, i) + MbI * Qins(iwlfi, 1) * sob_vsla + MfI * Qins(iwlfi, 1) * sof_vsla; 
        
        Qfhemi(iwlfo, i) = Qfhemi_v1(iwlfo, i) + Qfhemi_kt(iwlfo, i);
        
    else      % Multiple Scattering
        Qfdir(:, :, i) = Qsig(:, i) .* wleaf * rho;
        Qfdir(iwlfo, :, i) = Qfdir(iwlfo, :, i) + MQ * rho; 
        
        Qfhemi(:, i) = Qsig(:, i) .* wleaf * rho2;
        Qfhemi(iwlfo, i) = Qfhemi(iwlfo, i) + MQ * rho2; 
    end
    
    Qfyld(:, i) = MQ; 
    
    Qdown(:, i) = Qsig(:, i) .* wleaf * rho2;         
    Qdown(iwlfo, i) = Qdown(iwlfo, i) + MQ * rho2;
    
    % Prepare source for next collision order
    if i < 11
        Qsig(:, i + 1) = Qsig(:, i) .* wleaf * p;             
        Qsig(iwlfo, i + 1) = Qsig(iwlfo, i + 1) + MQ * p;
    end
end

% -------------------------------------------------------------------------
% 2. Black Soil (BS) Problem - Diffuse Illumination
% -------------------------------------------------------------------------
Qfdir_d = zeros(nb, angle, 11);    
Qfhemi_d = zeros(nb, 11); 
Qfdir_d_v1 = zeros(nb, angle, 11);
Qfdir_d_kt = zeros(nb, angle, 11);
Qfhemi_d_v1 = zeros(nb, 11);
Qfhemi_d_kt = zeros(nb, 11);

Qapar_d = zeros(nb, 11);   
Qdown_d = zeros(nb, 11);  
Qsig_d  = zeros(nb, 11);     
Qfyld_d = zeros(nf, 11);       

% First collision with incident diffuse light
Qsig_d(:, 1) = Qind * id;                

for i = 1:11
    Qapar_d(:, i) = Qsig_d(:, i) .* (1 - wleaf);  
    MQ = Mf * Qsig_d(iwlfi, i);
    
    if i == 1
        % Sunlit Crown
        Qfdir_d_v1(:, :, i) = Qind .* rho_l * sob_vsla_dif + Qind .* tau_l * sof_vsla_dif;    
        Qfdir_d_v1(iwlfo, :, i) = Qfdir_d_v1(iwlfo, :, i) + MbI * Qind(iwlfi, 1) * sob_vsla_dif + MfI * Qind(iwlfi, 1) * sof_vsla_dif; 
        Qfdir_d_v1(iwlfo, :, i) = Kc(iwlfo, :) .* Qfdir_d_v1(iwlfo, :, i);  
        
        % Shaded Crown
        Qfdir_d_kt(:, :, i) = Qind .* rho_l * sob_vsla_dif + Qind .* tau_l * sof_vsla_dif;    
        Qfdir_d_kt(iwlfo, :, i) = Qfdir_d_kt(iwlfo, :, i) + MbI * Qind(iwlfi, 1) * sob_vsla_dif + MfI * Qind(iwlfi, 1) * sof_vsla_dif; 
        Qfdir_d_kt(iwlfo, :, i) = Kt(iwlfo, :) .* sqrt(Ps_dir_inKz(iwlfo, :)) .* Qfdir_d_kt(iwlfo, :, i);  
        
        Qfdir_d(iwlfo, :, i) = Qfdir_d_v1(iwlfo, :, i) + Qfdir_d_kt(iwlfo, :, i);
        
        % Hemispherical components
        Qfhemi_d_v1(:, i) = Qind .* rho_l * sob_vsla_hemi_dif + Qind .* tau_l * sof_vsla_hemi_dif;
        Qfhemi_d_v1(iwlfo, i) = Qfhemi_d_v1(iwlfo, i) + MbI * Qind(iwlfi, 1) * sob_vsla_hemi_dif + MfI * Qind(iwlfi, 1) * sof_vsla_hemi_dif;    
        
        Qfhemi_d_kt(:, i) = Qind .* rho_l * sob_vsla_hemi_dif + Qind .* tau_l * sof_vsla_hemi_dif;
        Qfhemi_d_kt(iwlfo, i) = Qfhemi_d_kt(iwlfo, i) + MbI * Qind(iwlfi, 1) * sob_vsla_hemi_dif + MfI * Qind(iwlfi, 1) * sof_vsla_hemi_dif;    
        
        Qfhemi_d(iwlfo, i) = Qfhemi_d_v1(iwlfo, i) + Qfhemi_d_kt(iwlfo, i);
    else
        Qfdir_d(:, :, i) = Qsig_d(:, i) .* wleaf * rho;
        Qfdir_d(iwlfo, :, i) = Qfdir_d(iwlfo, :, i) + MQ * rho; 
        
        Qfhemi_d(:, i) = Qsig_d(:, i) .* wleaf * rho2;
        Qfhemi_d(iwlfo, i) = Qfhemi_d(iwlfo, i) + MQ * rho2; 
    end    
    
    Qfyld_d(:, i) = MQ; 
    
    Qdown_d(:, i) = Qsig_d(:, i) .* wleaf * rho2;         
    Qdown_d(iwlfo, i) = Qdown_d(iwlfo, i) + MQ * rho2;
    
    if i < 11
        Qsig_d(:, i + 1) = Qsig_d(:, i) .* wleaf * p;             
        Qsig_d(iwlfo, i + 1) = Qsig_d(iwlfo, i + 1) + MQ * p;
    end
end

% Sum up direct and diffuse radiation for the BS problem
Qapar_bs = sum(Qapar + Qapar_d, 2);
Qfdir_bs = sum(Qfdir + Qfdir_d, 3);   
Qfyld_bs = sum(Qfyld + Qfyld_d, 2);

% -------------------------------------------------------------------------
% 3. Preparation for the Soil (S) Problem
% -------------------------------------------------------------------------
Qdown_bs = Qins * t0 + Qind * td + sum(Qdown + Qdown_d, 2);   
Qind_s   = Qdown_bs .* rs;   

Qdown_bs_hot = Qins * t0;
Qind_s_hot   = Qdown_bs_hot .* rs;   

Qdown_bs_d = Qind * td + sum(Qdown + Qdown_d, 2);
Qind_s_d   = Qdown_bs_d .* rs;

% -------------------------------------------------------------------------
% 4. Soil (S) Problem - Soil-Canopy Multiple Interactions
% -------------------------------------------------------------------------
Qfdir_s  = zeros(nb, angle, 11); 
Qfhemi_s = zeros(nb, 11);   
Qapar_s  = zeros(nb, 11);   
Qdown_s  = zeros(nb, 11);  
Qsig_s   = zeros(nb, 11);     
Qfyld_s  = zeros(nf, 11);       

Qapar_ss = zeros(nb, 1);
Qfdir_ss = zeros(nb, angle);
Qfyld_ss = zeros(nf, 1);

for k = 1:8    % Loop for multiple bounces between soil and canopy
    
    if k == 1
        Qsig_s(:, 1) = Qind_s_hot * id + Qind_s_d * id;    
    else
        Qsig_s(:, 1) = Qind_s * id;     
    end
    
    for i = 1:11
        Qapar_s(:, i) = Qsig_s(:, i) .* (1 - wleaf);  
        
        MQ = Mf * Qsig_s(iwlfi, i);
        Qfdir_s(:, :, i) = Qsig_s(:, i) .* wleaf * rho;
        Qfdir_s(iwlfo, :, i) = Qfdir_s(iwlfo, :, i) + MQ * rho; 
        
        Qfhemi_s(:, i) = Qsig_s(:, i) .* wleaf * rho2;
        Qfhemi_s(iwlfo, i) = Qfhemi_s(iwlfo, i) + MQ * rho2; 
        
        Qfyld_s(:, i) = MQ; 
        
        Qdown_s(:, i) = Qsig_s(:, i) .* wleaf * rho2;         
        Qdown_s(iwlfo, i) = Qdown_s(iwlfo, i) + MQ * rho2;
        
        if i < 11
            Qsig_s(:, i + 1) = Qsig_s(:, i) .* wleaf * p;             
            Qsig_s(iwlfo, i + 1) = Qsig_s(iwlfo, i + 1) + MQ * p;
        end
    end
    
    Qapar_ss = Qapar_ss + sum(Qapar_s, 2);
    
    if (k == 1)
        Qfdir_ss = Qfdir_ss + sum(Qfdir_s, 3) + ...
                   Kc .* (Qins .* rs * kg + Qind .* rs * kg_dif) + ...
                   Kg .* (Qins + Qind) .* rs + ...
                   Kz .* (Qins + Qind) .* sqrt(Ps_dir_inKz) .* rs + ...
                   Kt .* ((Qins + Qind) * kg_kt .* rs) + ...
                   sum(Qdown + Qdown_d, 2) .* rs * tv;   
    else
        Qfdir_ss = Qfdir_ss + sum(Qfdir_s, 3) + Qind_s * tv;                                                                      
    end
    
    Qfyld_ss = Qfyld_ss + sum(Qfyld_s, 2);
    Qdown_ss = sum(Qdown_s, 2);
    
    Qind_s = Qdown_ss .* rs;  
end

% -------------------------------------------------------------------------
% 5. Final Output Assembly
% -------------------------------------------------------------------------
Qfdir_all = Qfdir_bs + Qfdir_ss;      
Qfyld_all = Qfyld_bs + Qfyld_ss;      
Qapar_all = Qapar_bs + Qapar_ss;      

end