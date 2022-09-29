function pardot = freeavg(t,par,m,paramp,Kinf,zt_lin)
% this is the ode function called by RK45_and_Avging.m, using the method of
% averaging for steady state response. This function calculates the rate of
% change of amplitude and phase, which are integrated to obtain the
% amplitude and phase.

% pulling out Iwan parameters from arrays
F_S = paramp(:,1);
Kt = paramp(:,2);             
chi = paramp(:,3);
beta = paramp(:,4);
    phi_max = F_S.*(1+beta)./(Kt.*(beta + (chi+1)./(chi+2)));
    R = F_S.*(chi+1)./((beta+ (chi+1)./(chi+2)).*phi_max.^(chi+2) );

    
% calculating rate of change of amplitude and phase
D = 4 * R/((chi + 3)*(chi + 2)) * par(1)^(chi + 3);                         %dissipation energy per cycle
r = par(1)/phi_max;
Kj = Kt * (1 - (r^(chi + 1)/(chi + 2)/(1 + beta)));                         %Fmax/X
wn = sqrt((Kj + Kinf)/m);                                                  %natural angular velocity(/frequency)
zeta = D/(m*2*pi*wn^2*par(1)^2)+zt_lin;                                     %damping ratio, taking wd = wn
pardot(1,1) = -wn *zeta*par(1);
pardot(2,1) = wn*sqrt(1 - zeta^2);

end