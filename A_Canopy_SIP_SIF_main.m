% =========================================================================
% Canopy-SIP SIF Model Demo
% 
% DESCRIPTION:
% This script demonstrates the simulation of canopy solar-induced 
% chlorophyll fluorescence (SIF) anisotropy. It combines the geometric-
% optical (GO) theory and spectral invariants theory (p-theory).
%
% AUTHOR:
% Yachang He, Yelu Zeng, Dalei Hao, ...
% Email: akhyc13@gmail.com ; zengyelu@cau.edu.cn
% Date: 2026-02-25
%
% REFERENCE:
% If you use this code in your research, please cite the following paper:
% He, Y., Zeng, Y., Hao, D., Shabanov, N. V., Huang, J., Yin, G., ... & Rossini, M. (2025). 
% "Combining geometric-optical and spectral invariants theories for modeling 
% canopy fluorescence anisotropy." Remote Sensing of Environment, 323, 114716.
% DOI: https://doi.org/10.1016/j.rse.2025.114716
% =========================================================================

clear all;
clc

%% 0. Initialize Wavelengths & Load Data
wls = (400:2561)';  % Optical wavelength range [nm]
wlf = (640:850)';   % SIF wavelength range [nm]
wle = (400:750)';   % PAR wavelength range [nm]

[~, iwlfi] = intersect(wls, wle);
[~, iwlfo] = intersect(wls, wlf);
[~, iwlfd] = intersect(wls, (710:790)');

nb = length(wls);
nf = length(iwlfo);
ne = length(iwlfi);

% Load input spectra and CI data
load('leafopt_P6.mat');         % Fluspect output （叶片模拟结果，输入参数）
load('Esun_SIP.mat');           % Solar irradiance
load('rs_01.mat');              % Soil reflectance

% Input the structural information of the scene (calculated by LESS model)
load('./CI_2/gap_tot.mat');     % Total gap fraction
load('./CI_2/gap_within.mat');  % Within-crown gap fraction
load('./CI_2/gap_betw.mat');    % Between-crown gap fraction
load('./CI_2/CI_within.mat');   % Within-crown clumping index
PATH_root = 'HET01_true_all.txt'; % Used for Hemispherical Interceptance (iD) calculation

%% 1. Prepare for Simulation
% 1.1 Sun-Sensor Geometry
angle = 25;  % Number of simulation angles
va = zeros(angle, 4); % View angle matrix [VZA, VAA, -, -]

% Set up specific view angles (e.g., along the principal plane)
for t = 1:12
    va(t, 1) = 5 * (13 - t); % View Zenith Angle (VZA)
    va(t, 2) = 0;            % View Azimuth Angle (VAA) - Forward
end
% va(13, :) remains 0, representing the Nadir view
for t = 1:12
    va(t + 13, 1) = 5 * t;   % View Zenith Angle (VZA)
    va(t + 13, 2) = 180;     % View Azimuth Angle (VAA) - Backward
end

SZA = 20;  % Sun Zenith Angle [degrees]
SAA = 0;   % Sun Azimuth Angle [degrees]

% 1.2 Vegetation Type (Discrete  Canopy)

% Gap fraction from CI (Index 13 represents Nadir, 9 represents Solar angle)
gap_H        = gap_betw(13, 3);      % Nadir between-crown gap fraction
gap_H_within = gap_within(13, 3);    % Nadir within-crown gap fraction
gap_H_tot    = gap_tot(13, 3);       % Nadir total gap fraction
CI_H_within  = CI_within(13, 3);     % Nadir within-crown clumping index

gap_S        = gap_betw(9, 3);       % Solar direction between-crown gap fraction
gap_S_within = gap_within(9, 3);     % Solar direction within-crown gap fraction
gap_S_tot    = gap_tot(9, 3);        % Solar direction total gap fraction


% Scene information for HET01
Crowndeepth = 12.86;  % Average crown depth [m] 
Height      = 20;       % Canopy height [m]
Height_c    = 6.87;   % Crown center height [m]
dthr        = 0.86;  % Diameter to Height Ratio
bl          = 0.2;      % Leaf width [m]
HotSpotPar  = 0.02;     % Hotspot parameter at leaf scale

c1 = CI_H_within;       % Within-crown clumping index (Nadir)
c2 = 1 - gap_H; 
c  = c1 * c2;           % Canopy-level clumping index

iD   = getHemiInterceptancev4_H(PATH_root); % Scene hemispherical interceptance
LAI  = 5;               % Single crown sphere LAI [m2/m2]
FAVD = 0.375;           % Foliage Area Volume Density
D    = 0;               % Ratio of diffuse to incoming irradiance

% 1.3 Soil and Canopy Properties (LIDF setup)
TypeLidf = 2; % 1 = Two-parameter description, 2 = Single-parameter (Campbell)

if (TypeLidf == 1)    
    % Two-parameter LIDF: |LIDFa| + |LIDFb| < 1
    LIDFa = 0;
    LIDFb = 0;
    [lidf, litab] = dladgen(LIDFa, LIDFb);   
elseif (TypeLidf == 2)
    % Single-parameter LIDF (Average leaf angle in degrees)
    % Planophile (26.76), Erectophile (63.24), Spherical (57.3), etc.
    LIDFa = 57.3;   
    LIDFb = 0;
    [lidf, litab] = campbell(LIDFa);        
end
lidf = lidf';

% Leaf and soil optics
rho = zeros(2162, 1);
tau = zeros(2162, 1);
rho(1:2101, 1) = leafopt.refl'; % Leaf reflectance
tau(1:2101, 1) = leafopt.tran'; % Leaf transmittance 
w  = rho + tau;                 % Leaf single scattering albedo
rg = rs;                        % Soil reflectance  

% Structural parameters for modeling
Znum = floor(Crowndeepth ./ HotSpotPar .* FAVD); % Number of layers
go_par = dthr * Crowndeepth;                     % Hotspot parameter at crown scale

%% 2. Start Simulation
% -------------------------------------------------------------------------
MfI  = leafopt.MfI;
MfII = leafopt.MfII;
MbI  = leafopt.MbI;
MbII = leafopt.MbII;
    
% Incident radiation
Qins = result_SIP(:, 2);        % Direct solar irradiance 
Qins = Qins * cosd(SZA);        % Adjust for sun zenith angle
Qind = Qins .* 0;               % Diffuse irradiance (set to 0)

% adjustment factor for Photosystem I
MfI = MfI * 1.29;
MbI = MbI * 1.29;

% Core Radiative Transfer and SIF Calculations
[data_p_rho] = get_p_rho_comment(LAI, SZA, SAA, va, angle, lidf, Height, Height_c, Crowndeepth, gap_tot, gap_betw, gap_within, gap_S, gap_S_tot, gap_S_within, go_par, c, HotSpotPar, iD, D);

[Qfdir_1, Qfyld_1, Qapar_all] = get_rta_withsif_comment(nb, nf, ne, iwlfi, iwlfo, angle, Qins, Qind, rho, tau, rs, MfI, MbI, data_p_rho);
[Qfdir_2, Qfyld_2, Qapar_all] = get_rta_withsif_comment(nb, nf, ne, iwlfi, iwlfo, angle, Qins, Qind, rho, tau, rs, MfII, MbII, data_p_rho);
[Qpdir_all, Qpyld_all, Qpapar_all, Qpdir_bs] = get_rta_nosif_comment(nb, nf, ne, iwlfi, iwlfo, angle, Qins, Qind, rho, tau, rs, data_p_rho);

% Calculate SIF directional emission by subtracting background
SRTE_Fs_fdir1 = Qfdir_1 - Qpdir_all;
SRTE_Fs_fdir1 = SRTE_Fs_fdir1(iwlfo, :)./pi;

SRTE_Fs_fdir2 = Qfdir_2 - Qpdir_all;
SRTE_Fs_fdir2 = SRTE_Fs_fdir2(iwlfo, :)./pi;

SRTE_Fs_fdir_all = SRTE_Fs_fdir1 + SRTE_Fs_fdir2;  % Total SIF emission
SRTE_RefAll = Qpdir_all ./ repmat((Qins + Qind), 1, angle);

%% 3. Save Results
% Save the simulation results for GitHub demo
save('demo_result_2026.mat', 'SRTE_Fs_fdir_all', 'SRTE_Fs_fdir1', 'SRTE_Fs_fdir2');
disp('Simulation completed successfully. Results saved to demo_result_2026.mat');


%% 4. Plot Results
% 绘制天底方向 (Nadir view, 对应 index 13) 的 SIF 光谱图
disp('Generating SIF spectrum plot...');

figure('Name', 'Canopy SIF Spectrum', 'Color', 'w');
nadir_idx = 13; % 第13个角度为天底方向 (VZA = 0)

% 绘制总 SIF、PS I 和 PS II
plot(wlf, SRTE_Fs_fdir_all(:, nadir_idx), 'r-', 'LineWidth', 2);


% 图表格式化 (满足学术出版基本要求)
xlabel('Wavelength (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('SIF Radiance (W m^{-2} sr^{-1} um^{-1})', 'FontSize', 12, 'FontWeight', 'bold'); 
title('Simulated Canopy SIF Spectrum at Nadir View', 'FontSize', 14);
legend('Total SIF', 'Location', 'northeast');
xlim([640 850]);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1);
box on;
hold off;

disp('Plot generated successfully.');