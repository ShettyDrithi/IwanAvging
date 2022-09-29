function xdot = forcingRK(t_RK45,x_RK45,m,ffunc,paramp,zt_lin,Kinf,funclogic)
% this is the ode function called by RK45_and_Avging.m when there is a
% forcing term in the EoM.

persistent zt wn X; 
% NOTE: persistent variables have been defined to use the previous value
% every time the function is called by ode45.

% pulling out Iwan parameters from arrays
F_S = paramp(1);
Kt = paramp(2);             chi = paramp(3);
beta = paramp(4);
    phi_max = F_S.*(1+beta)./(Kt.*(beta + (chi+1)./(chi+2)));
    R = F_S.*(chi+1)./((beta+ (chi+1)./(chi+2)).*phi_max.^(chi+2) );


% initializing the persistent variables the first time the function is called
if isempty(zt) || t_RK45 == 0
    zt = zt_lin;
end
if isempty(wn) || t_RK45 == 0
    wn = sqrt((Kt + Kinf)/m);
end
if isempty(X) || t_RK45 == 0
    X = max(eps,sqrt(x_RK45(1)^2 + (x_RK45(2)/wn)^2));
end

if funclogic
    fext = ffunc(t_RK45);
else
    ts_pos = min([(floor(t_RK45/ffunc{1})+1),length(ffunc{2})-1]);
    fext = ffunc{2}(ts_pos)*(1-rem(t_RK45/ffunc{1},1)) + ffunc{2}(ts_pos+1).*rem(t_RK45/ffunc{1},1);
%     fext = interp1(ffunc{1},ffunc{2},t_RK45);
end
dif = 1;
count = 1;
count_max = 100;

while (abs(dif)>1e-13) && count < count_max
    if count ~= 1
        D = 4 * R/((chi + 3)*(chi + 2)) * X^(chi + 3);                     % dissipation energy per cycle
        r = X/phi_max;
        Kj = Kt * (1 - (r^(chi + 1)/(chi + 2)/(1 + beta)));                % Fmax/X
        wn = sqrt((Kj + Kinf)/m);                                          % natural angular velocity(/frequency)
        zt = D/(m*2*pi*wn^2*X^2)+zt_lin;
        X = sqrt(x_RK45(1)^2 + x_RK45(2)^2/(wn*sqrt(1-zt^2))^2);
    end
    Xprev = X;
    D = 4 * R/((chi + 3)*(chi + 2)) * X^(chi + 3);                         
    r = X/phi_max;
    Kj = Kt * (1 - (r^(chi + 1)/(chi + 2)/(1 + beta)));                
    wn = sqrt((Kj + Kinf)/m); 
    zt = D/(m*2*pi*wn^2*X^2)+zt_lin;
    X = sqrt(x_RK45(1)^2 + x_RK45(2)^2/(wn*sqrt(1-zt^2))^2);
    dif = X - Xprev;
    count = count + 1;
    if count == count_max
        error('RK45 max. iteration count reached')
    end
end

xdot(1,1) = x_RK45(2);
xdot(2,1) = fext/m - 2*zt*wn*x_RK45(2) - wn^2*x_RK45(1);
end