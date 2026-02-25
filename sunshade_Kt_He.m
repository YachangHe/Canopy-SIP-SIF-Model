function [kc, kg] = sunshade_Kt_He(tts, tto, psi, Gs, Go, CIs, CIo, LAI)
% SUNSHADE_KT_HE Calculates sunlit foliage and soil probabilities (No Hotspot).
%
% DESCRIPTION:
% This function computes the bi-directional gap probabilities assuming spatial 
% independence between the solar and viewing directions (i.e., Pso = Ps * Po). 
% It does NOT include the hotspot correlation enhancement. This specific 
% formulation is typically utilized for the shaded crown component (Kt) in 
% the continuous canopy scenario.
%
% INPUTS:
%   tts - Sun zenith angle [degrees]
%   tto - View zenith angle [degrees]
%   psi - Relative azimuth angle [degrees]
%   Gs  - G-function value in the solar direction
%   Go  - G-function value in the viewing direction
%   CIs - Clumping index in the solar direction
%   CIo - Clumping index in the viewing direction
%   LAI - Leaf Area Index [m2/m2]
%
% OUTPUTS:
%   kc  - Probability of viewing sunlit foliage (independent assumption)
%   kg  - Probability of viewing sunlit soil (independent assumption)

% -------------------------------------------------------------------------
deg2rad = pi / 180;

cos_tts = cos(tts * deg2rad);             
tan_tto = tan(tto * deg2rad);             
cos_tto = cos(tto * deg2rad);             
tan_tts = tan(tts * deg2rad);             

% Ensure relative azimuth angle is symmetric around the principal plane
psi = abs(psi - 360 * round(psi / 360));  

% -------------------------------------------------------------------------
% Case 1: Exact overlap direction
if tts == tto && psi == 0
    kc = 1 - exp(-Gs * CIs / cos_tts * LAI);
    kg = exp(-Gs * CIs / cos_tts * LAI);
else
% Case 2: General bi-directional calculation without hotspot correlation
    
    nl = 20; % Number of canopy layers for numerical integration
    x  = (-1/nl : -1/nl : -1)';   % Column vector representing layer depths
    xl = [0; x];                  % Add top of canopy level
    dx = 1/nl;
    iLAI = LAI / nl;              % LAI of each elementary layer
    
    k = Gs / cos_tts; % Optical path length factor (Solar)
    K = Go / cos_tto; % Optical path length factor (View)
    
    % Independent probabilities of viewing a leaf
    Ps = exp(k * xl * CIs * LAI);                                              
    Po = exp(K * xl * CIo * LAI);                                              
    
    % Correct Ps/Po for finite layer thickness (dx)
    Ps(1:nl) = Ps(1:nl) .* (1 - exp(-k * CIs * LAI * dx)) / (k * CIs * LAI * dx);                                      
    Po(1:nl) = Po(1:nl) .* (1 - exp(-K * CIo * LAI * dx)) / (K * CIo * LAI * dx);  
    
    Pso = zeros(size(Po));
    
    % Numerically integrate the joint independent probability
    for j = 1:length(xl)
        Pso(j, :) = quad(@(y) Psofunction(K, k, CIs, CIo, LAI, y), xl(j) - dx, xl(j)) / dx; %#ok<FREMO>
    end
    
    % Take care of rounding errors
    Pso(Pso > Po) = min([Po(Pso > Po), Ps(Pso > Po)], [], 2);    
    Pso(Pso > Ps) = min([Po(Pso > Ps), Ps(Pso > Ps)], [], 2);    
    
    % Calculate final visible foliage and soil
    kc = iLAI * CIo * K * sum(Pso(1:nl));  
    kg = Pso(nl + 1);                      
end

end

% =========================================================================
% INTERNAL FUNCTION: Independent joint probability (No hotspot)
% =========================================================================
function pso = Psofunction(K, k, CIs, CIo, LAI, xl)
    % Calculates joint probability assuming independence: Pso = Ps * Po
    pso = exp((K * CIo + k * CIs) * LAI * xl); 
end