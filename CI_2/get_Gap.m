% 本函数的目的是根据 CI 求 方向间隙率

function Gap = get_Gap(VZA,CI,LAI,G)

theta = deg2rad(VZA);

Gap = exp(-CI * G * LAI / cos(theta));

end