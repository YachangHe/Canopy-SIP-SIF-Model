function iD = getHemiInterceptancev4_H(PATH_root)

% This routine calculates hemeispherical interceptance
% 计算半球拦截概率
% by Yachang
% 对每个方向的总间隙率 积分就行

%PATH_root = 'HET01_true_all.txt';
Gap_tot = readmatrix(PATH_root, 'NumHeaderLines', 1);    % 总间隙率

%%
za = [0 : 5 : 85 89];

iv = ones(size(za));

for i = 1:size(za,2)

        p = Gap_tot(i,3);        
        iv(i) = (1-p) * sind(2.*za(i));
end

iD = trapz(deg2rad(za), iv);    % required to be checked!!!

end

