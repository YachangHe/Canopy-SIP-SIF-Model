function g_L = get_gFun( iorien,theta_L)
%This routine calculates leaf angle distribution gL(thetaL)
%type of leaf normal orientation:
%1-planophile (a=1,b=2); 2-erectophile (a=-1,b=2);
%3-plagiophile (a=-1, b=4); 4-extremophile(a=1,b=4);
%5-uniform (a=0, b any);6-spherical(g=sin(theta))
%theta_L in radian

%%
if(iorien==6)     %6-spherical
    g_L=sin(theta_L);
else
    switch(iorien)
        case 1     %1-planophile
            a=1;b=2;
        case 2     %2-erectophile
            a=-1;b=-2;
        case 3      %3-plagiophile
            a=-1;b=4;
        case 4      %4-extremophile
            a=1;b=4;
        case 5      % 5-uniform
            a=0;b=0;
    end
    g_L=(2/pi)*(1+a*cos(b*theta_L));    
end

