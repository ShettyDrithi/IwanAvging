%% example script to run iwanavg.m 
clear;clc;
% close all;
%% Build the Model
%
% /|                           |--> u1
% /|                       _________
% /|         Kinf         |         |
% /|______/\/\/\/\/\______|         |
% /|   ________________   |   M0    |--->f_ext
% /|__| Fs Kt beta chi |__|         |
% /|  |________________|  |_________|
% /|

%% Define the system parameters:

m = 1;                  % Mass
Kinf = 1.4212e+05;      % linear stiffness
zt_material = 1e-4;     % modal damping ratio
wninf = sqrt(Kinf/m);   % in rad/s
sysprop = [m,Kinf,zt_material];

% Nonlinear Iwan Joint parameters
Fs   = 400;          % slip force
% Fs   = 40000;          % slip force
Kt  = 2.5266e+05;      % stiffness of the joint
beta  = 1;             % parameter that is associated with energy dissipation curve
chi = -0.5;            % parameter associate with slope of log(energy dissipation) vs. log(force)
paramp = [Fs,Kt,chi,beta];
wn0 = sqrt((Kinf+Kt)/m);% frequency when joint fully stuck in rad/s
fn0 = wn0/2/pi;

% Force parameters
% T_forcePulse= 0.020;
% T_forcePulse= 1/fn0;
% Am = 200;            % Impact Force Levels
% tsim = 20;             % simulation duration, seconds
% Force = [T_forcePulse Am];
% ffunc=@(t)Am*sin(pi/T_forcePulse*t).*(t<T_forcePulse); %defining impulse as a func of time
% Tband = [0.1 15];

% properties for sine beat (uncomment to try sine beat)
F = 100;                % freq for max amplitude
dF = 10; 
tsim = linspace(0,10,20000);
t0 = 0;
FAmp = 100;
ffunc = @(time)ping_create(F,dF,t0,FAmp,time);
Tband = [0.4,10];
figure; plot(tsim,ffunc(tsim));title('force vs time')

[resp_avg,simtime_avg,var_avg] = iwanavg(ffunc,tsim,paramp,m,Kinf,zt_material,'Tband',Tband);
% pulling out the required variables from the varargout structure
    x_avg = resp_avg(1,:); 
    v_avg = resp_avg(2,:); 
    wn_avg = var_avg.wn_fit;
    zt_avg = var_avg.zt_fit;
    vfit_avg = var_avg.yfit;
    vfit_avg = abs(vfit_avg);
    Y_avg = var_avg.Yv;
    w_avg = var_avg.w;
    t_avg = resp_avg(4,:);
   
%% calculating the damping and frequency when fext ~=0 using amplitude approx. algorithm
% {
% find the last time instant having nonzero force
fext = ffunc(t_avg);
t_endind = find(abs(fext)>0,1,'last'); 
xforce = x_avg(1:t_endind);
vforce = v_avg(1:t_endind);

    phi_max = Fs.*(1+beta)./(Kt.*(beta + (chi+1)./(chi+2)));
    R = Fs.*(chi+1)./((beta+ (chi+1)./(chi+2)).*phi_max.^(chi+2) );
    
% initializing ---- these are the vectors you want!!!!!!!
Xdforce = zeros(t_endind,1);
wnforce = zeros(t_endind,1);
Xforce = zeros(t_endind,1);
wdforce = zeros(t_endind,1);
ztforce = zeros(t_endind,1);

X = sqrt(xforce(1)^2 + vforce(1)^2/(wn0^2*(1 - zt_material^2)));
for j = 1:t_endind
    if X == 0
        wn = wn0;
        zt = zt_material;
    else
        dif = 1;
        count = 1;
        count_max = 200;
        while (abs(dif)>1e-13) && count < count_max
            D = 4 * R/((chi + 3)*(chi + 2)) * X^(chi + 3);
            r = X/phi_max;
            Kj = Kt * (1 - (r^(chi + 1)/(chi + 2)/(1 + beta)));
            wn = sqrt((Kj + Kinf)/m);
            zt = D/(m*2*pi*wn^2*X^2)+zt_material;
            wd = wn*sqrt(1-zt^2);
            Xprev = X;
            X = sqrt((xforce(j))^2 + vforce(j)^2/(wd)^2);
            dif = X - Xprev;
            count = count + 1;
            if count == count_max
                error('max. iteration count reached')
            end
        end
    end
Xforce(j) = X;
Xdforce(j) = X*wd;
wdforce(j) = wd;
wnforce(j) = wn;
ztforce(j) = zt;
end
%}

% plots
figure(24)
semilogx(Xdforce,ztforce,'-b')
xlabel('|velocity|');ylabel('damping \zeta')
title('obtained using closed-form expressions no Hilbert')

figure(25)
semilogx(Xdforce,wnforce,'-b')
xlabel('|velocity|');ylabel('frequency \omega_n')
title('obtained using closed-form expressions no Hilbert')
    
%% plotting response
figure(11)
plot(t_avg,x_avg,'-r')
ylabel('response displacement');xlabel('time')

% plotting amplitude dependent damping and natural frequency
figure(12)
semilogx(vfit_avg,zt_avg,'-r')
xlabel('response amplitude');ylabel('damping \zeta')

figure(13)
semilogx(vfit_avg,wn_avg,'-r')
xlabel('response amplitude');ylabel('natural frequency \omega_n')

% plotting damping and frequency obtained directly from iwanavg algo vs
% time
figure(14)
semilogy(tsim(1:t_endind),ztforce,var_avg.tfree,var_avg.ztnohilb,'-b')
xlabel('time');ylabel('damping \zeta')
title('obtained using closed-form expressions no Hilbert')

figure(15)
semilogy(tsim(1:t_endind),wnforce,var_avg.tfree,var_avg.wnnohilb,'-b')
xlabel('time');ylabel('natural frequency \omega_n')
title('obtained using closed-form expressions no Hilbert')
                                                            