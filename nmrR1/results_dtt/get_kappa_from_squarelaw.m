%-------------------------------------------------------------------------
% get_kappa_from_squarelaw.m
%
% Function that calculates the bending rigidity from the slope of the
% square law dependence
%
% -------------------------------------------------------------------------
% Input parameters
%
% data      two-column vector with Scd^2 and R1
% mpp       full bilayer thickness in Angstrom
% lfeq      Larmor frequency
% T         temperature in Kelvin
% plotflag  flag to plot the data and fit (1) or not (0)  
%
% -------------------------------------------------------------------------
% Output parameters
%
% kappa     bending rigidity in units of kT
% -------------------------------------------------------------------------

function [kappa] = get_kappa_from_squarelaw(data, mpp, lfeq, T, plotflag)

eta = 0.1; % bilayer viscosity in Pascal * second
Ss = 0.5; % slow order parameter, unitless
t = mpp*(10^(-10)); % thickness in meters
kT = 1.380649*(10^(-23))*T; % kbT in J

% find slope of square-law dependence
[f, g] = fit(data(:,1),data(:,2),'poly1');
f = coeffvalues(f);
slope = f(1); % must be in s^(1/2)

chiQ = 167*(10^-3)*(10^6); % in 1/s
omega0 = lfeq*(10^6); % in 1/s
D = slope / ( (15/8)*(pi^2)*(chiQ^2)*sqrt(1/omega0) );

if plotflag==1
    figure
    hold on
    plot(data(:,1),data(:,2),'o');

    x=0:0.001:max(data(:,1));
    y=f(1)*x+f(2);
    plot(x,y);
end

K = ((9*(kT^2)*eta)./(50*(pi^2)*(Ss.^4)*(D^2))).^(1/3); % in Newton

kappa = (K*t)./kT; % kT
