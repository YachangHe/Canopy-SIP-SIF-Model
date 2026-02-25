% 用来测试解析的CI和间隙率; 

Gs = 0.6;  % G函数
VZA = 0;
R = 10 ;   % 球体半径

leaf_s = 1000;   % 总叶面积
LAI = leaf_s / (pi*R^2);      % 叶面积指数
LAD = leaf_s / (4*pi*R^3/3);   % 叶面积体密度

% 定义积分函数
f = @(x) x .* exp(-LAD * Gs * x); 
% 计算定积分从 0 到 20
result = integral(f, 0, 2*R);

disp(['总叶面积: ', num2str(leaf_s)]);
% 显示结果
Gap_result = result./(2*R^2);     % 间隙率的解析计算结果
disp(['Gap解析结果: ', num2str(Gap_result)]);
theta = deg2rad(VZA);
CI = ( -cos(theta) * log(Gap_result)) / (Gs * LAI) ;   % CI的解析计算结果
disp(['CI解析结果: ', num2str(CI)]);