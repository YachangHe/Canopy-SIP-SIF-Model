%% 这个程序用来测试天底方向聚集指数


Gap_tot = 0.97229;  %  单颗树场景总间隙，by Less

Gap_betw = 1 - (pi * 10^2 ) / (100^2) ;   % 单颗树场景大间隙，可以直接计算

Gap_within = Gap_tot - Gap_betw;  % 冠层内的间隙率 (相对于大场景)

Gap_within = (Gap_within*100*100)./(pi*10^2);   % 冠层内的间隙率 (相对于冠层球体)

LAI_within = 5 ;
iorien = 5;        % 叶倾角类型
%type of leaf normal orientation(iorien):
%1-planophile (a=1,b=2); 2-erectophile (a=-1,b=2);
%3-plagiophile (a=-1, b=4); 4-extremophile(a=1,b=4);
%5-uniform (a=0, b any);6-spherical(g=sin(theta))

% num_rows = size(Gap_tot, 1);
% CI_tot = zeros(num_rows, 1);
% CI_within = zeros(num_rows, 1);
%for i = 1:num_rows

    TypeLidf = 1;
    if (TypeLidf==1)    % 双参数描述，不要取的太极限，0.99表示1
    %	LIDF type 		a 		 b
    %	Planophile 		1		 0
    %	Erectophile    -1	 	 0
    %	Plagiophile 	0		-1
    %	Extremophile 	0		 1
    %	Spherical 	   -0.35 	-0.15
    %	Uniform         0        0
    % 	requirement: |LIDFa| + |LIDFb| < 1	
    LIDFa	=	0;
    LIDFb	=	0;
   
    elseif (TypeLidf==2)
    %	Planophile 		26.76
    %	Erectophile     63.24
    %	Plagiophile 	45
    %	Extremophile 	45
    %	Spherical 	    57.3
    %	Uniform         45
    LIDFa	=	45;   % average leaf angle (degrees) 0 = planophile	/	90 = erectophile
    LIDFb	=	0;
    end

    if (TypeLidf==1)
    [lidf,litab] = dladgen(LIDFa,LIDFb);   % 双参数描述
    elseif (TypeLidf==2)
    [lidf,litab] = campbell(LIDFa);        % 单参数描述
    end
    lidf=lidf';


    VZA = 0;
    VAA = 0;
    
    [Gss,kk]    =   PHASE2(VZA,lidf);   % by sail 计算的G函数
    %LAI_within = 15 * 500 * pi / (100*100*(1 - Gap_betw(i,3)));  %树冠内LAI=单叶面积 /((1 - 刚体球间隙率)* 场景总面积)  方向LAI
    % 这个地方需要再确认下，CI分母的LAI是带方向的吗？
    %LAI_within = 5;  % 不带方向性的LAI，单个球体的LAI
    Gv = get_G(iorien, VZA, VAA);  % 聚集指数的G，与观测方向有关   by 伟华 计算的G函数
    %CI_tot= get_CI(VZA, Gap_tot(i, 3), LAI_crown, Gv);
    
%     Gv = 0.66;
%     tao1 = Gv*0.375;
%     Gap_within_true = (1-(2*tao1*10+1)^(-2*tao1*10))/(2*tao1^2*10^2);
%     CI_within_true = -log(Gap_within_true) / (Gv * 5);
%     运行 66 -69 lines 发现  CI_within_true 始终无法计算得到 2/3。 这跟谁的G函数更准已经没有关系了。
%     本质上 比拼的是谁的间隙率算的更准，因为别的都是一样的。   2024.09.13

    CI_within_PATH = get_CI(VZA, Gap_within, LAI_within, Gv);
    CI_within_sail = get_CI(VZA, Gap_within, LAI_within, Gss);
% end
