% 用来测试解析的CI和间隙率; 

% Gs = 0.6;  % G函数
% VZA = 0;
% R = 10 ;   % 球体半径
% 
% leaf_s = 2500;   % 总叶面积
% LAI = leaf_s / (pi*R^2);      % 叶面积指数
% LAD = leaf_s / (4*pi*R^3/3);   % 叶面积体密度
% 
% % 定义积分函数
% f = @(x) x .* exp(-LAD * Gs * x); 
% % 计算定积分从 0 到 20
% result = integral(f, 0, 2*R);
% 
% disp(['LAI: ', num2str(LAI)]);
% % 显示结果
% Gap_result = result./(2*R^2);     % 间隙率的解析计算结果
% disp(['Gap解析结果: ', num2str(Gap_result)]);
% theta = deg2rad(VZA);
% CI = ( -cos(theta) * log(Gap_result)) / (Gs * LAI) ;   % CI的解析计算结果
% disp(['CI解析结果: ', num2str(CI)]);



% % 测试Less 
% LAD = 0.3;
% Gs = 0.6;  % G函数
% VZA = 0;
% R = 10 ;   % 球体半径
% LAI = LAD * (4*pi*R^3/3) / (pi*R^2);
% Gap_tot = 0.97405;
% Gap_betw = 1 - (pi * 10^2 ) / (100^2) ;   % 单颗树场景大间隙，可以直接计算
% Gap_within = Gap_tot - Gap_betw;  % 冠层内的间隙率 (相对于大场景)
% Gap_result = (Gap_within*100*100)./(pi*10^2);   % 冠层内的间隙率 (相对于冠层球体)
% 
% theta = deg2rad(VZA);
% CI = ( -cos(theta) * log(Gap_result)) / (Gs * LAI) ;   % CI的解析计算结果
% disp(['CI解析结果: ', num2str(CI)]);
% % CI_gap_test_byYachang
% % CI解析结果: 0.72866
% 
% 
% % 测试Less 
% LAD = 0.6;
% Gs = 0.6;  % G函数
% VZA = 0;
% R = 10 ;   % 球体半径
% LAI = LAD * (4*pi*R^3/3) / (pi*R^2);
% Gap_tot = 0.97027;
% Gap_betw = 1 - (pi * 10^2 ) / (100^2) ;   % 单颗树场景大间隙，可以直接计算
% Gap_within = Gap_tot - Gap_betw;  % 冠层内的间隙率 (相对于大场景)
% Gap_result = (Gap_within*100*100)./(pi*10^2);   % 冠层内的间隙率 (相对于冠层球体)
% 
% theta = deg2rad(VZA);
% CI = ( -cos(theta) * log(Gap_result)) / (Gs * LAI) ;   % CI的解析计算结果
% disp(['CI解析结果: ', num2str(CI)]);
% % CI_gap_test_byYachang
% % CI解析结果: 0.60937


% 测试Less 
LAD = 0.6;
Gs = 0.5;  % G函数
VZA = 0;
R = 10 ;   % 球体半径
LAI = LAD * (4*pi*R^3/3) / (pi*R^2);
Gap_tot = 0.97027;
Gap_betw = 1 - (pi * 10^2 ) / (100^2) ;   % 单颗树场景大间隙，可以直接计算
Gap_within = Gap_tot - Gap_betw;  % 冠层内的间隙率 (相对于大场景)
Gap_result = (Gap_within*100*100)./(pi*10^2);   % 冠层内的间隙率 (相对于冠层球体)
theta = deg2rad(VZA);
CI = ( -cos(theta) * log(Gap_result)) / (Gs * LAI) ;   % CI的解析计算结果
disp(['CI模拟结果: ', num2str(CI)]);
% CI_gap_test_byYachang
% CI解析结果: 0.80374



% 定义积分函数
f = @(x) x .* exp(-LAD * Gs * x); 
% 计算定积分从 0 到 20
result = integral(f, 0, 2*R);

% 显示结果
Gap_result = result./(2*R^2);     % 间隙率的解析计算结果
disp(['Gap解析结果: ', num2str(Gap_result)]);
theta = deg2rad(VZA);
CI = ( -cos(theta) * log(Gap_result)) / (Gs * LAI) ;   % CI的解析计算结果
disp(['CI解析结果: ', num2str(CI)]);