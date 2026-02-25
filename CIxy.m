function [CIs]=CIxy(CIy1,CIy2,tts)

% CIXY Calculates the angular dependence of the clumping index.
% DESCRIPTION:
% This function estimates the clumping index at a specific zenith angle 
% using a linear interpolation model. It accounts for the phenomenon that 
% the canopy gap fraction and clumping effect change with the viewing or 
% illumination angle.
%
% INPUTS:
%   CIy1 - Clumping index at nadir (zenith angle = 0 degrees)
%   CIy2 - Clumping index at a large zenith angle (typically 75 degrees)
%   tts  - Target zenith angle (solar or view zenith angle) [degrees]
%
% OUTPUTS:
%   CIs  - Estimated clumping index at the target zenith angle (tts)

% if tts<20
%     CIs=CIy1;
% else if tts>60
%     CIs=CIy2;   
%     else
% CIs=(CIy2-CIy1)/(60-20)*(tts-20)+CIy1;
%     end
% end


% -------------------------------------------------------------------------
% Empirical linear interpolation
CIs = (CIy2 - CIy1) / (75 - 0) * (tts - 0) + CIy1;

% Physical boundary check: Clumping index should not increase indefinitely
% at extreme zenith angles. Cap it at CIy2 (or 1.0) for tts > 75.
if tts > 75
    CIs = CIy2; 
end

end