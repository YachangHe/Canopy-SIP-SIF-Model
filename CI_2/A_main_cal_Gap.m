%% 这个程序用来计算方向间隙率，根据输入的CI 计算 gap
%  类比lwh的模型，他用路径长度得到gap；本研究用CI得到gap

% function Gap = get_Gap(VZA,CI,LAI,G)
% 读取.mat文件 
% 这个code 用来 读取聚集指数，并计算间隙率
clear all
clc

filename_1 = 'CI_tot.mat';
load(filename_1);    % 总聚集指数

filename_2 = 'CI_within.mat';
load(filename_2);   % 冠层内聚集指数


LAI_crown = 2.356;          % 已知LAI; 场景总的LAI
iorien = 5;        % 叶倾角类型
%type of leaf normal orientation(iorien):
%1-planophile (a=1,b=2); 2-erectophile (a=-1,b=2);
%3-plagiophile (a=-1, b=4); 4-extremophile(a=1,b=4);
%5-uniform (a=0, b any);6-spherical(g=sin(theta))

num_rows = size(CI_tot, 1);
gap_tot = zeros(num_rows, 1);
gap_within = zeros(num_rows, 1);
gap_within_2 = zeros(num_rows, 1);
gap_betw = zeros(num_rows, 1);

for i = 1:num_rows
    VZA = CI_tot(i, 1);
    VAA = CI_tot(i, 2);
    %LAI_within = 15 * 500 * pi / (100*100*(1 - Gap_tot(i,3)));  %树冠内LAI=单叶面积 /((1 - 刚体球间隙率)* 场景总面积) 
    LAI_within = 5;
    Gv = get_G(iorien, VZA, VAA);  % 聚集指数的G，与观测方向有关

    % get_Gap(VZA,CI,LAI,G) 
    gap_tot(i) = get_Gap(VZA, CI_tot(i, 3), LAI_crown, Gv);    % 总间隙率
    gap_within(i) = get_Gap(VZA, CI_within(i,3), LAI_within, Gv);  %  注意这个地方算的是冠层内的间隙率 （投影面积对应的是球体投影面积，不是全场景）
    gap_within_2(i) = gap_within(i).*(pi*100*15)./100./100;          %  冠层内间隙率.*投影面积 ./ 场景总面积
    gap_betw(i) = gap_tot(i) - gap_within_2(i);   % 大间隙
end

% 拼接数据
gap_tot = [CI_tot(:, 1:2), gap_tot];
gap_within = [CI_tot(:, 1:2), gap_within];
gap_betw = [CI_tot(:, 1:2), gap_betw];
% 保存结果到新的.mat文件
save('gap_tot.mat', 'gap_tot');
save('gap_within.mat', 'gap_within');
save('gap_betw.mat', 'gap_betw');
% 保存到新的.mat文件
%save('processed_data.mat', 'abs_third_col');
