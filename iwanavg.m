function [resp,simtime,varargout] = iwanavg(Forcing,time,paramp,m,Kinf,zt_lin,varargin)
% This function computes the response of a SDOF system with 1 Iwan joint, 
% in the micro-slip region, using the closed form expressions for 
% dissipation and stiffness derived by Segalman. The Hilbert transform is 
% then used to obtain the natural frequency and damping ratio.
%
% /|                           |--> u1
% /|                       _________
% /|         Kinf         |         |
% /|______/\/\/\/\/\______|         |
% /|   ________________   |   m     |--->f_ext
% /|__| Fs Kt beta chi |__|         |
% /|  |________________|  |_________|
% /|
% [algorithm explained in D. Shetty, M. Allen, ASME JVA 2020]

% Input parameters:
% Forcing : externally applied force as a function of time
%           or externally applied force vector. In this case, time
%           vetor should have constant spacing.
%           If simulating free response, pass in initial amplitude through varargin          
% time : Simulation time vector (if a vector isn't provided, an equally
%        spaced vector having dt=1/(20*fn0) will be created)
% paramp : physical Iwan parameters [Fs,Kt,chi,beta]
% m : mass
% Kinf,zt_lin : stiffness of elastic spring, linear damping ratio
%
% **optional varargin: enter name-value pairs (NOT case-sensitive)
% Tband : time band for the piecewise-linear hilbert analysis, default []
% ICs : initial displacement and velocity entered as an array, default [0,0]
%
% Output parameters:
% resp : [displacement;velocity;accn;time]
%            matrix with the firt row being the response displacement,
%            second row being the velocity, third acceleration and fourth, time
% simtime : time taken for numerical integration
%
% **optional varargout:
% wn_fit,zt_fit,yfit : natural frequency, damping ratio and response 
%                      velocity amplitude obtained from Hilbert transform 
% y,w : complex FFT coefficients of response velocity, frequency vector
%       corresponding to 'y' in rad/s
% 
% note: this function calls subfunctions forcingRK.m and freeavg.m to seperately solve the transient and 
%       free response problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parsing varargin for additional input
pObj = inputParser;
addParameter(pObj,'TBand',[]);
addParameter(pObj,'ICs',[]);

parse(pObj,varargin{:});
TBand = pObj.Results.TBand;
ICs = pObj.Results.ICs;


% pulling out the Iwan parameters from the array and converting them to
% their mathematical equivalents
Fs = paramp(1);
Kt = paramp(2);
chi = paramp(3);
beta = paramp(4);

fn0 = sqrt(Kt+Kinf)/2/pi;           % stuck natural frequency

% creating a time vector if only simulation period given
if length(time) == 1
    dt = 1/(fn0*20);                % 20 times arbitrarily chosen and found sufficient
    time = (0:dt:time);  
end

% ensure time is a row vector
if size(time,1)>1
    time = time.';
end

% assigning the forcing function or intial amplitude
if isa(Forcing,'function_handle')
    ffunc = Forcing;
    fext=ffunc(time);
    funclogic = 1;

else
    fext = Forcing;
    dt = time(2) - time(1);
    ffunc = {dt,fext};
end

% initial conditions for ode
if ~isempty(ICs)
    init = ICs;
else
    init = [eps,0];
    % NOTE: the initial amplitude should ideally be 0, but the equations
    % lead to a divide by zero error and hence fail.
end
% finding the index of the last non-zero element in the force vector
t_endind = find(abs(fext)>0,1,'last');  
if isempty(t_endind)
    t_endind = 1;
end

if t_endind ~= length(time)
    Nstp_rk = find(time>=1/fn0,1,'first');
else
    Nstp_rk = 0;
end

tic
    tspan_rk = time(1:(t_endind+Nstp_rk)); 
%     tspan_rk = time(1:(t_endind+round(0.05*length(time)))); 
%   NOTE: 1% of time instants for which the force is present is added to
%   ensure entire transient response is integrated. 1% here is arbitrarily
%   chosen, scope for improvement
    init_rk = init;
%     options_rk = odeset('RelTol',1e-10,'AbsTol',1e-10);
%   NOTE: observed here that having a tolerance when external force is
%   present does not significantly improve the accuracy but results in
%   major computation burden when the force is present for longer periods
%   of time. Therefore using default tolerance values for forced response

    [t_rk,x_rk] = ode45(@(t_rk,x_rk)forcingRK(t_rk,x_rk,m,ffunc,...
                      paramp,zt_lin,Kinf,funclogic),tspan_rk,init_rk);
t_rk = t_rk.';
x_rk = x_rk.';
a_rk = gradient(x_rk(2,:),t_rk(2) - t_rk(1));

%     sol_rk = ode45(@(t_rk,x_rk)forcingRK(t_rk,x_rk,m,ffunc,...
%                       paramp,paramm,zt_lin,Kinf,A,funclogic),tspan_rk,init_rk);
% t_rk = tspan_rk;
% [x_rk,xdot_rk] = deval(sol_rk,tspan_rk);
% a_rk = (xdot_rk(2,:));

if length(tspan_rk) < length(time)
% calculating the amplitude and phase at the end of ode45 using multsineobj
% {    
options = optimset('TolFun',1e-10,'TolX',1e-10);    
    fn_guess = fn0;
    T_begin = 1/fn_guess;
    n_step = ceil(T_begin/(t_rk(2)-t_rk(1)));
fopt = fminsearch(@(x) multsineobj(x,x_rk(1,(end-n_step):end),...
                    t_rk((end-n_step):end)),fn_guess,options);
    T_opt = 1/fopt;
    n_step = ceil(T_opt/(t_rk(2)-t_rk(1)));
    % overwriting the above optimization if the result is beyond the length
    % of the time vector since in that case, the T_begin value is more
    % accurate than the optimized result
    % (this needs some improvement)
    if n_step>=length(t_rk)
        n_step = ceil(T_begin/(t_rk(2)-t_rk(1)));
    end
    [~,~,cs_amp] = multsineobj(fopt,x_rk(1,(end-n_step):end),t_rk((end-n_step):end));
    
    comp_amp=[cs_amp(1); cs_amp(2:2:end)-1i*cs_amp(3:2:end)];
        
        phase_rk = angle(comp_amp(2)*exp(1i*2*pi*fopt*t_rk(end)));
        X_rk = abs(comp_amp(2));
%}


tspan_avg = time((t_endind+Nstp_rk):end);
% tspan_avg = time((t_endind+round(0.05*length(time))):end);
init_avg = [X_rk phase_rk];
options_avg = odeset('RelTol',1e-12,'AbsTol',1e-12);
%   NOTE: for the free response, using a low tolerance value to get higher sccuracy 
%   when response diminishes to nearly zero values. This gives better
%   amplitude-dependent damping and frequency results.
Teval_avg = time((t_endind+Nstp_rk+1):end);
% Teval_avg = time((t_endind+round(0.05*length(time))+1):end);

sol = ode45(@(t,par)freeavg(t,par,m,paramp,Kinf,zt_lin),...
                                        tspan_avg,init_avg,options_avg);
    % t here is a dummy variable - nothing is assigned to it.
simtime = toc;

[par,pardot] = deval(sol,Teval_avg);

% finding the NL damping and frequency (needed for accn)
phi_max = Fs.*(1+beta)./(Kt.*(beta + (chi+1)./(chi+2)));
R = Fs.*(chi+1)./((beta+ (chi+1)./(chi+2)).*phi_max.^(chi+2) );

D = 4 * R/((chi + 3)*(chi + 2)) .* real(par(1,:)).^(chi + 3);                         
r = real(par(1,:))/phi_max;
Kj = Kt * (1 - (r.^(chi + 1)/(chi + 2)/(1 + beta)));                         
wn = sqrt((Kj + Kinf)/m); 
zeta = D./(m*2*pi*wn.^2.*real(par(1,:)).^2)+zt_lin;  
% zeta = sqrt((D./(m*2*pi*wn.^2.*real(par(1,:)).^2)).^2+zt_lin^2) + zt_lin;  

x_avg = real(par(1,:) .* exp(par(2,:)*1i));
v_avg = real(pardot(1,:).*exp(par(2,:)*1i) + par(1,:).*exp(par(2,:)*1i).*(pardot(2,:)*1i));
a_avg = real(-2*zeta.*wn.*v_avg - wn.^2.*x_avg);
    xt = [x_rk(1,:), x_avg];
    vt = [x_rk(2,:), v_avg];
    at = [a_rk, a_avg];
    ts = [t_rk, Teval_avg];
else
    simtime = toc;
    xt = x_rk(1,:);
    vt = x_rk(2,:);
    at = a_rk;
    ts = t_rk;
end
resp = [xt; vt; at; ts];

[Yv,ws] = fftp(vt,(ts(2)-ts(1)));
[Ya,~] = fftp(at,(ts(2)-ts(1)));

% Hilbert transform only when Tband is provided
    if ~isempty(TBand)
%         [wn_fit,zt_fit,~,yfit]=hilbssm_edited(vt,ts,20,[],Tband); 
        [wn_fit,zt_fit,~,yfit]=hilblsmMP(vt,ts,21,'setBand',TBand,'noPlots');
    else 
        wn_fit = [];zt_fit = []; yfit = [];
    end
    
if nargout > 2
    var.wn_fit = wn_fit;
    var.zt_fit = zt_fit;
    var.yfit = yfit;
    var.Yv = Yv;
    var.Ya = Ya;
    var.w = ws;
    var.wnnohilb = wn;
    var.ztnohilb = zeta;
    var.tfree = Teval_avg;
    var.vfree = v_avg;
    varargout{1} = var;
end