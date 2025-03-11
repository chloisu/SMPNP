phi_wall = 100

phi_centerline = 2;
phi = @(x) phi_wall*exp(-10*x)+phi_centerline;

a0 = 1;
aa = 0.5;ca
ac = 0.5;
ca_centerline = 1;
cc_centerline = 2;
xa_centerline = ca_centerline*aa^3;
xc_centerline = cc_centerline*ac^3;
Ca =(aa/a0)^3.*log(xa_centerline/(1-xa_centerline-xc_centerline)) - phi_centerline;
Cc =(ac/a0)^3.*log(xc_centerline/(1-xa_centerline-xc_centerline)) + phi_centerline;

Ac = @(x) exp((Cc-phi(x))*(a0/ac)^3);
Aa = @(x) exp((Ca+phi(x))*(a0/aa)^3);

ca = @(x) 1.0/aa^3.*(Aa(x)-Aa(x).*Ac(x)./(1+Ac(x)))./(1+Aa(x)-Aa(x).*Ac(x)./(1+Ac(x)));
cc = @(x) 1.0/ac^3.*(Ac(x)-Aa(x).*Ac(x)./(1+Aa(x)))./(1+Ac(x)-Aa(x).*Ac(x)./(1+Aa(x)));

x = linspace(0,1);
%plot(x,phi(x))
hold on
%plot(x,Ac(x))
plot(x,ca(x))
plot(x,cc(x))