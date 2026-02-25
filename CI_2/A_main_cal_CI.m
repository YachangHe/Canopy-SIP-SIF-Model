%% 输入间隙率，这个程序用来计算方向聚集指数
% function CI = get_CI(VZA,Gap,LAI,G)
% 读取.txt文件 从第二行开始读取
% 这个code 用来 读取 less 计算的间隙率，并计算聚集指数

filename_1 = 'Gap_true_05.txt';
Gap_tot = readmatrix(filename_1);    % 总间隙率

filename_2 = 'Gap_ball_05.txt';
Gap_betw = readmatrix(filename_2);   % 大间隙率

Gap_within = Gap_tot(:,3)-Gap_betw(:,3);                  %     冠层内的间隙率  by Less  相对于大场景

Gap_within = (Gap_within*100*100)./(pi*100*15);          %  冠层内的间隙率  by Less  相对于单个球体

LAI_crown = 2.356;          % 已知LAI; 场景总的LAI   Rami 单个球体的LAI = 5
iorien = 5;        %  叶倾角类型

%type of leaf normal orientation(iorien):
%1-planophile (a=1,b=2); 2-erectophile (a=-1,b=2);
%3-plagiophile (a=-1, b=4); 4-extremophile(a=1,b=4);
%5-uniform (a=0, b any);6-spherical(g=sin(theta))

num_rows = size(Gap_tot, 1);
CI_tot = zeros(num_rows, 1);
CI_within = zeros(num_rows, 1);

for i = 1:num_rows
    VZA = Gap_tot(i, 1);
    VAA = Gap_tot(i, 2);
    %LAI_within = 15 * 500 * pi / (100*100*(1 - Gap_betw(i,3)));  %树冠内LAI=单叶面积 /((1 - 刚体球间隙率)* 场景总面积)  方向LAI
    LAI_within = 5;  % 不带方向性的LAI，单个球体的LAI
    Gv = get_G(iorien, VZA, VAA);  % 聚集指数的G，与观测方向有关
    CI_tot(i) = get_CI(VZA, Gap_tot(i, 3), LAI_crown, Gv);
    CI_within(i) = get_CI(VZA, Gap_within(i), LAI_within, Gv);
end

% 拼接数据
CI_tot = [Gap_tot(:, 1:2), CI_tot];
CI_within = [Gap_tot(:, 1:2), CI_within];

% 保存结果到新的.mat文件
save('CI_tot.mat', 'CI_tot');
save('CI_within.mat', 'CI_within');

% 保存到新的.mat文件
%save('processed_data.mat', 'abs_third_col');
