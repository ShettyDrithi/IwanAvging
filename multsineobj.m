function [value,x_fit,cs] = multsineobj(f,x,t,varargin)
%  Represent a measured signal as a Fourier Series and return the
%  coefficients as well as a measure of how well it fits the measurement.
%
%  [residual,y_fit,cs] = multsineobj(f,y,t,N);
%
%  input:       f   - fundamental frequency
%               y   - measured signal
%               t   - time sequence
%               varargin  - number of harmonics
%  output:      value     - sum of squared deviations
%               y_fit     - fourier series fit results
%               cs        - amplitudes (coefficients) of the following set
%                       of functions, in this order:
%       [DC, cos(1*omega*t), sin(1*omega*t), cos(2*omega*t), sin(2*omega*t),
%       ...cos(N*omega*t),sin(N*omega*t)]
%
% This routine is meant to be used as the % objective function to an 
% optimization routine to find the precise fundamental frequency of a 
% signal that is a sum of sinusoids that are integer multiples of the 
% fundamental.
%
%   For example, to fit a multi-sine with 4 frequencies to the data y:
%
%     options = optimset('TolFun',1e-8,'TolX',1e-8);
%     f0 = 10; % frequency of initial guess, in Hz
%     fopt = fminsearch(@(x) multsineobj(x,y,t,4),f0,options)
%
%   Finds the frequency assuming that y is a sum of sines and
%   cosines out to the 4th harmonic on the time vector t.
%
% Then, to see a plot:
%
%     [value,y_fit,cs] = multsineobj(fopt,y,t,4);
%     figure(2); plot(t,y,t,y_fit);
%
% To get a vector of complex amplitudes:
%   X=[cs(1); cs(2:2:end)-1i*cs(3:2:end)];
%
% M.S. Allen, July 2009
%
%

if length(x)~=numel(x); error('x must be a vector'); end
if length(t)~=numel(t); error('t must be a vector'); end
x = x(:);
t = t(:);

if nargin < 4
A = [ones(size(t)), cos(f*2*pi*t), sin(f*2*pi*t),cos(2*f*2*pi*t), sin(2*f*2*pi*t),...
    cos(3*f*2*pi*t), sin(3*f*2*pi*t),cos(4*f*2*pi*t), sin(4*f*2*pi*t)];
else
    N = varargin{1};
    if length(N) == 1;
        A = zeros(length(t),1+2*N);
        A(:,1) = 1;
        for k = 1:N
            A(:,2*k:2*k+1) = [cos(k*f*2*pi*t), sin(k*f*2*pi*t)];
        end
    else
        % This signifies that a vector of harmonics is given.
        % Fit a multi-sine only at those frequencies.
        ns = N;
        N = length(ns);
        A = zeros(length(t),1+2*N);
        A(:,1) = 1;
        for k = 1:N
            A(:,2*k:2*k+1) = [cos(ns(k)*f*2*pi*t), sin(ns(k)*f*2*pi*t)];
        end
    end
end
b = x;

cs = A\b;
x_fit = A*cs;

value = sum((b - x_fit).^2);

