function [kc, kg] = sunshade_H(tts, tto, psi, Gs, Go, CIs, CIo, LAI, hotspot)
% SUNSHADE_H Calculates sunlit foliage and soil probabilities (Direct Hotspot Input).
%
% DESCRIPTION:
% This function computes the bi-directional gap probabilities considering the 
% hotspot effect based on the Kuusk (1995) formulation. 
% 
% NOTE: Unlike the standard 'sunshade' function which estimates the hotspot 
% parameter from leaf size and canopy height, this version takes the hotspot 
% parameter directly as an input. This approach is consistent with the PATH 
% radiative transfer model.
%
% INPUTS:
%   tts     - Sun zenith angle [degrees]
%   tto     - View zenith angle [degrees]
%   psi     - Relative azimuth angle [degrees]
%   Gs      - G-function value in the solar direction
%   Go      - G-function value in the viewing direction
%   CIs     - Clumping index in the solar direction
%   CIo     - Clumping index in the viewing direction
%   LAI     - Leaf Area Index [m2/m2]
%   hotspot - Direct hotspot parameter (e.g., related to crown/leaf scale)
%
% OUTPUTS:
%   kc      - Probability of viewing sunlit foliage (sunlit crown fraction)
%   kg      - Probability of viewing sunlit soil (bi-directional gap fraction)

% -------------------------------------------------------------------------
deg2rad = pi / 180;

cos_tts = cos(tts * deg2rad);             
tan_tto = tan(tto * deg2rad);             
cos_tto = cos(tto * deg2rad);             
tan_tts = tan(tts * deg2rad);             

% Ensure relative azimuth angle is symmetric around the principal plane
psi = abs(psi - 360 * round(psi / 360));  

% -------------------------------------------------------------------------
% Case 1: Exact hotspot direction (Sun and View directions overlap perfectly)
if tts == tto && psi == 0
    kc = 1 - exp(-Gs * CIs / cos_tts * LAI);
    kg = exp(-Gs * CIs / cos_tts * LAI);
else
% Case 2: General bi-directional calculation using canopy stratification
    
    nl = 20; % Number of canopy layers for numerical integration
    x  = (-1/nl : -1/nl : -1)';   % Column vector representing layer depths
    xl = [0; x];                  % Add top of canopy level
    dx = 1/nl;
    iLAI = LAI / nl;              % LAI of each elementary layer
    
    % Use directly inputted hotspot parameter
    q = hotspot;
    
    % Angular distance between solar and viewing directions
    dso = sqrt(tan_tts.^2 + tan_tto.^2 - 2 * tan_tts .* tan_tto .* cos(psi * deg2rad));
    
    k = Gs / cos_tts; % Optical path length factor (Solar)
    K = Go / cos_tto; % Optical path length factor (View)
    
    % Independent probabilities of viewing a leaf
    Ps = exp(k * xl * CIs * LAI);                                              
    Po = exp(K * xl * CIo * LAI);                                              
    
    % Correct Ps/Po for finite layer thickness (dx)
    Ps(1:nl) = Ps(1:nl) .* (1 - exp(-k * CIs * LAI * dx)) / (k * CIs * LAI * dx);                                      
    Po(1:nl) = Po(1:nl) .* (1 - exp(-K * CIo * LAI * dx)) / (K * CIo * LAI * dx);  
    
    Pso = zeros(size(Po));
    
    % Numerically integrate the joint probability over each canopy layer
    for j = 1:length(xl)
        Pso(j, :) = quad(@(y) Psofunction(K, k, CIs, CIo, LAI, q, dso, y), xl(j) - dx, xl(j)) / dx; %#ok<FREMO>
    end
    
    % Take care of rounding errors
    Pso(Pso > Po) = min([Po(Pso > Po), Ps(Pso > Po)], [], 2);    
    Pso(Pso > Ps) = min([Po(Pso > Ps), Ps(Pso > Ps)], [], 2);    
    
    % Calculate final visible sunlit foliage and sunlit soil
    kc = iLAI * CIo * K * sum(Pso(1:nl));  % Visible sunlit leaf
    kg = Pso(nl + 1);                      % Visible sunlit soil
end

end

% =========================================================================
% INTERNAL FUNCTION: Kuusk hotspot correlation function
% =========================================================================
function pso = Psofunction(K, k, CIs, CIo, LAI, q, dso, xl)
    if dso ~= 0
        alf = (dso / q) * 2 / (k + K);
        pso = exp((K * CIo + k * CIs) * LAI * xl + sqrt(K * CIo * k * CIs) * LAI / alf * (1 - exp(xl * alf)));
    else
        pso = exp((K * CIo + k * CIs) * LAI * xl - sqrt(K * CIo * k * CIs) * LAI * xl); 
    end
end