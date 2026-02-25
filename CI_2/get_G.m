function G_Fun = get_G( iorien,theta,fi )
%This routine estimates integral G(Omega)
%(1/2*pi)int(0,pi/2)int(0,2*pi)g_L |(Omega,Omega_L)|d theta_L d fi
%Here g_L(Omega_L)=(2/pi)(1+a*cos(b*theta_L)

%type of leaf normal orientation(iorien):
%1-planophile (a=1,b=2); 2-erectophile (a=-1,b=2);
%3-plagiophile (a=-1, b=4); 4-extremophile(a=1,b=4);
%5-uniform (a=0, b any);6-spherical(g=sin(theta))
%theta,fi indicate the direction Omega in Degree !

%%
theta = deg2rad(theta);
fi = deg2rad(fi);

n=30;
m=4*n;
h_theta=0.5*pi/n;
h_fi=2*pi/m;
theta_i=0.5*h_theta;
fi_1=0.5*h_fi;
G_Fun=0;
for i=1:n
    fi_j=fi_1;
    c_i=cos(theta_i);
    s_i=sin(theta_i);
    xx=0;
    for j=1:m
        yy=cos(theta)*c_i+sin(theta)*s_i*cos(fi-fi_j);
        xx=xx+abs(yy);
        fi_j=fi_j+h_fi;
    end  
    xx=xx*h_fi;
    yy=get_gFun(iorien,theta_i);
    G_Fun=G_Fun+yy*xx;
    theta_i=theta_i+h_theta;
end
G_Fun=h_theta*G_Fun;
G_Fun=G_Fun/(2*pi);
G_Fun_cos = G_Fun./cos(theta);

end

