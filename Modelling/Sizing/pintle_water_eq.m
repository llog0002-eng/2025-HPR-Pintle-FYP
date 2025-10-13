clear all; close all; clc;

% 3.5kN Liquid Engine sizing script for N2O/IPA
% Written by Philip Pine
% Updated by Luke Logan
% Last Updated 16/04/2025

% References:  https://ir.library.oregonstate.edu/concern/defaults/bv73c785v?locale=en
% Thermodynamic performance from NASA CEA & Rocket Propulsion Analysis
% (RPA)
% Fuel density from REFPROP 

% Constants
go = 9.81; % gravitational acceleration (m/s^2)

% Design points
T = 3000; % thrust (N)
Pc = 25; % desired chamber pressure (bar)
Ps_f = 12.65; % fuel supply pressure (bar)
Ps_ox = Ps_f; % oxidiser supply pressure (bar)
Pa = 1.01325; % sea level pressure (bar) 
OFtarg = 2.5; % desired oxidiser to fuel ratio

% Values from NASA CAE at given OF and chamber pressure for N2O/IPA at design pressure 

AeAt = 3.7627; % exit area to nozzle throat area from NASA CEA

Ce = 2039.4; % effective exhaust velocity (m.s) from NASA CEA
Cstar = 1433.7; % characteristic velocity (m/s) (propellant combustion property) from NASA CEA
isp = Ce/go; % specific impulse (s)

% Properties of propellant
% Temperature of fluids at injector
Tinj = 278; % [K]

% Density of water (replacing liquid N2O)
rho_ox = 998; % [kg/m^3]
% Density of water (replacing gaseous N2O)
rho_ox_g = 998; % [kg/m^3]
% Viscosity of water (replacing liquid N2O)
mu_ox = 1.002e-3; % [Poise]
% Surface tension of water (replacing liquid N2O)
sigma_ox = 72.8e-3; % [N/m]

% Density of water (replacing IPA)
rho_f = 998; % [kg/m^3]
% Viscosity of water (replacing IPA)
mu_f = mu_ox; % [Poise]
% Surface tension of water (replacing IPA)
sigma_f = sigma_ox; % [N/m]

gamma_t = 1.2772; % spefici heat ratio at the throat
gamma_g = 1.2739; % specific heat ratio of the exhaust 

cp_t = 1.9007; % specific heat at constant pressure at the throat kJ/kgK
cp_g = 1.9829; % specific heat at constant presssure at the exit of nozzle kJ/kgK

cv_t = cp_t/gamma_t;
cv_g = cp_g/gamma_g;

R_t = cp_t-cv_t; % gas constant at throat kJ/kgK
R_g = cp_g-cv_g; % gas cosntant for exhaust kJ/kgK

% Calculating required mass flow rates 
mp = (T/Ce); % required propellant mass flow rates (kg/s) % ASSUMES IDEALLY EXPANDED (SEA LEVEL OPERATION)
mftarg = (mp/(1+OFtarg)); % fuel flow rate (kg/s)
moxtarg = OFtarg*mftarg; % oxidiser mass flow rate (kg/s)

%% Nozzle sizing

At = Cstar*mp(end)/(Pc(end)*10^5); % nozzle throat area (m^2)
Dt = 2*(At/pi)^0.5; % diameter of throat (m)

Ae = AeAt*At; % nozzle exit area (m^2)
De = 2*(Ae/pi)^0.5; % diameter of exit (m)

theta = 15; % nozzle cone half angle (deg)
Ln = (De-Dt)/(2*tand(theta)); % nozzle length (m)

%% Chamber sizing 

% Calculating chamber volume from characterstic chamber length and throat area
Lstar = 1; % characteristic chamber length

Vch = Lstar*At; % chamber volume (m^3)

% Empirical formula for chamber area to throat area
et = 9; % (Ac/At) contraction ratio (ratio of cross sectional area of chamber divided by throat), from https://wikis.mit.edu/confluence/pages/viewpage.action?pageId=153816550
Ac = At*et; % calculating chamber area from throat area and contraction rato (m^2)
rc = (Ac/pi)^0.5; % radius of combustion chamber (m)
rt = Dt/2; % radius of throat (m)
Dc = 2*(Ac/pi)^0.5; % chamber diameter (m) %

Lc = Vch/Ac; % calculating combustion chamber length (m)

% conical contraction
theta_c = 45; 
Lfrus = (rc-rt)/(tand(theta_c)); % contraction length  (m)

Ltotal = Lfrus+Lc+Ln; % total length from combustion chamber to exit of nozzle (m)

%% Chamber/Nozzle Geometry Plotting

graphic = false;

if graphic == true
    f1 = figure();
    hold on
    grid on
    plot([0,Lc]*1e3,[Dc/2,Dc/2]*1e3,"k"); % Plot combustion chamber
    plot([Lc,Lc+Lfrus]*1e3,[Dc/2,Dt/2]*1e3,"k"); % Plot converging section
    plot([Lc+Lfrus,Lc+Lfrus+Ln]*1e3,[Dt/2,De/2]*1e3,"k"); % Plot nozzle
    daspect([1 1 1])
    xlabel("Length (mm)")
    ylabel("Radius (mm)")
    title("Combustion chamber + nozzle profile")
    f1.Position(3:4) = [1200,1200*De/Ltotal];
end

%% Pintle sizing
% Parameters
% DR - diameter ratio (chamber diameter and pintle diameter)
% SR - skip ratio
% Lopen - pintle opening distance

% Conditions for optimisation of pintle
% TMR close to 1
% Vaporisation length 
% Lopen > 0.10 m
% Pintle tip angle <= 20
% Blockage factor between 0.5 and 0.7 
% SR = 1 for maximum combustion efficiency

% TO MATCH:
% We = 150
% TMR = 1 (mox = 1.14, mf = 0.46)
% Re_o = 4860, Re_i = 8990

% Non dimensional pintle paramters
SR = 1.0; % skip distance ratio 
DR = 4.89202; % ratio of chamber diameter and pintle diameter (between 3-5) https://ltu.diva-portal.org/smash/get/diva2:1845405/FULLTEXT01.pdf

% Pintle geometry parameters 
% Dpt = Dc/DR*1e3; % pintle tip diameter (mm) - function of chamber diameter
Dpt = 25; % Override pintle diameter (mm)
Ls = SR*Dpt; % skip length (mm) - distance outer flow travels before impingement point

% Pressure difference over injector
dP_ox = (Ps_ox-Pa)*1e5; % (Pa)
dP_f = (Ps_f-Pa)*1e5; % (Pa)

% Throttle
throttle = 0.6; % 1 = full throttle

%% Annular gap pintle

t_sleeve = 5.5; % thickness of sleeve (mm)
id_sleeve = Dpt - 2*t_sleeve; % sleeve ID (mm)

% Discharge coefficients for inner and outer flows, from experimental data https://www.researchgate.net/publication/301440576_Experiments_with_Pintle_Injector_Design_and_Development
Cd_i = 0.7; % Inner orifice Cd
Cd_o = 0.7; % Outer orifice Cd
Cd_ip = 0.7; % Inner passthrough Cd
Cd_op = 0.7; % Outer passthrough Cd

theta_pt = 40; % pintle tip angle (deg, from horizontal)
Dpr = 3; % pintle rod diameter (mm) ## Change this to be a dependent variable later
Dcg = 4.5; % center gap diameter (mm) #s# Change this to be a dependent variable later
r_post = Dpt/2; % post diameter radius (mm)

%%% Correct pressures for passthrough holes
pass_in_d = 2.5; % Inner passthrough hole diameter (mm)
pass_in_n = 10; % Number of inner passthrough holes
A_ip = pass_in_n * pass_in_d^2/4 * pi; % Area of inner passthrough holes (mm2)

pass_o_d = 1.5; % Inner passthrough hole diameter (mm)
pass_o_n = 8; % Numberof inner passthrough holes
A_op = pass_o_n * pass_o_d^2/4 * pi; % Area of inner passthrough holes (mm2)

A_o = 18.6; % Outer orifice area (mm2)
A_i = 65.7*throttle; % Inner orifice area (mm2)

% Flow conductances
K_i = Cd_i * A_i / 1e3^2;
K_ip = Cd_ip * A_ip / 1e3^2;
K_o = Cd_o * A_o / 1e3^2;
K_op = Cd_op * A_op / 1e3^2;

K_i_eq = (1/K_i^2+1/K_ip^2)^(-1/2); % Flow conductance, combine in series
K_o_eq = (1/K_o^2+1/K_op^2)^(-1/2); % Flow conductance, combine in series
K_eq = K_i_eq + K_o_eq; % Flow conductance, combine in parallel

mox = K_i_eq * sqrt(2*rho_ox*Ps_ox*1e5);
mf = K_o_eq * sqrt(2*rho_f*Ps_f*1e5);
mp_pred = mox+mf;

OF = mox/mf; % Estimated OF ratio
Gap_o = sqrt(A_o/pi+(Dpt/2)^2)-Dpt/2; % Outer flow opening distance (mm)
Gap_i = A_i/(pi*id_sleeve); % Pintle opening distance (mm)
Gap_iz = Gap_i/cosd(theta_pt); % Pintle axial opening distance (mm)

Dh_i = 2*Gap_i/1000; % Hydraulic diameter for inner flow (m)
Dh_o = 2*Gap_o/1000; % Hydraulic diameter for outer flow (m)

U_i = mox/rho_ox/(A_i/(1e3)^2); % Velocity of inner flow (m/s)
U_o = mf/rho_f/(A_o/(1e3)^2); % Velocity of outer flow (m/s)

Cd_eff_i = mox/((A_i)*(1e-3)^2*sqrt(2*rho_f*Ps_f*1e5));
Cd_eff_o = mf/((A_o)*(1e-3)^2*sqrt(2*rho_f*Ps_f*1e5));
Cd_eff = (mox+mf)/((A_o+A_i)*(1e-3)^2*sqrt(2*rho_f*Ps_f*1e5));

%% Non-Dimensional Output Parameters

% Reynolds number
Re_i = rho_ox*U_i*Dh_i/mu_ox;
Re_o = rho_f*U_o*Dh_o/mu_f;

% Total momentum ratio
% Thanks chat for the derivation https://chatgpt.com/share/68ecc095-0eac-8013-bbe3-d4a45ccc6e4e
TMR = (mox*U_i)/(mf*U_o);

% Momentum flux ratio
J = (rho_f*U_o^2)/(rho_ox*U_i^2);

% Inner and outer Weber numbers
We_i = rho_ox_g*U_i^2*Gap_i*1e-3/sigma_ox;
We_o = rho_f*U_o^2*Gap_o*1e-3/sigma_f;
% Different inner and outer We determines spray morphology
% Want We_i or We_o ~ 3000 for "fully developed fan spray"
% "The spray and atomization process and its effects on combustion performance of pintle injector,” Ph.D. thesis

%% Spray Angle Prediction

% Half angle
theta_spray_half = atand(TMR * cosd(theta_pt)/(TMR * sind(theta_pt) + 1)); % Degrees
