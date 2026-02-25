function CI = get_CI(VZA,Gap,LAI,G)

theta = deg2rad(VZA);

CI = ( -cos(theta) * log(Gap)) / (G * LAI) ;

end

