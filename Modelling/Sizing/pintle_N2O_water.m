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
Ps_f = 5; % fuel supply pressure (bar)
Ps_ox = 50; % oxidiser supply pressure (bar)
Pa = 1.01325; % sea level pressure (bar) 
OF = 2.5; % desired oxidiser to fuel ratio

% Values from NASA CAE at given OF and chamber pressure for N2O/IPA at design pressure 

AeAt = 3.7627; % exit area to nozzle throat area from NASA CEA

Ce = 2039.4; % effective exhaust velocity (m.s) from NASA CEA
Cstar = 1433.7; % characteristic velocity (m/s) (propellant combustion property) from NASA CEA
isp = Ce/go; % specific impulse (s)

% Properties of propellant
% Temperature of fluids at injector
Tinj = 278; % [K]

% Density of liquid N2O from Coolprop @ T = 298 K and P = 34 bar
rho_ox = 1220; % [kg/m^3]
% Density of gaseous N2O from Coolprop @ T = 298 K and P = 34 bar
rho_ox_g = 77.0; % [kg/m^3]
% Viscosity of liquid N2O from wikipedia
mu_ox = 3.237e-3; % [Pa.s]
% Surface tension of liquid N2O at -25 C
sigma_ox = 0.0101; % [N/m]
% Gas constant of gaseous N2O
R_ox_g = 180; % [J/kg K]
% Saturation pressure of N2O
Pv_ox = Ps_ox; %

% Density of water
rho_f = 998; % [kg/m^3]
% Viscosity of water at 20 C
mu_f = 1.002e-3; % [Pa.s]
% Surface tension of water at room temp
sigma_f = 72.8e-3; % [N/m]

gamma_t = 1.2772; % speficic heat ratio at the throat
gamma_g = 1.2739; % specific heat ratio of the exhaust 

cp_t = 1.9007; % specific heat at constant pressure at the throat kJ/kgK
cp_g = 1.9829; % specific heat at constant presssure at the exit of nozzle kJ/kgK

cv_t = cp_t/gamma_t;
cv_g = cp_g/gamma_g;

R_t = cp_t-cv_t; % gas constant at throat kJ/kgK
R_g = cp_g-cv_g; % gas cosntant for exhaust kJ/kgK

% Calculating required mass flow rates 
mp = (T/Ce); % required propellant mass flow rates (kg/s) % ASSUMES IDEALLY EXPANDED (SEA LEVEL OPERATION)
mf_targ = (mp/(1+OF)); % fuel flow rate (kg/s)
mox_targ = OF*mf_targ; % oxidiser mass flow rate (kg/s)

%% Nozzle sizing

At = Cstar*mp(end)/(Pc(end)*10^5); % nozzle throat area (m^2)
Dt = 2*(At/pi)^0.5; % diameter of throat (m)

Ae = AeAt*At; % nozzle exit area (m^2)
De = 2*(Ae/pi)^0.5; % diameter of exit (m)

theta = 15; % nozzle cone half angle (deg)
Ln = (De-Dt)/(2*tand(theta)); % nozzle length (m)

%% Chamber sizing 

% Calculating chamber volume from characterstic chamber length and throat
% area
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
% Injection velocity (10-50 m/s) ?? sus

% Non dimensional pintle paramters
SR = 1.0; % skip distance ratio 
DR = 4.89202; % ratio of chamber diameter and pintle diameter (between 3-5) https://ltu.diva-portal.org/smash/get/diva2:1845405/FULLTEXT01.pdf

% Pintle geometry parameters 
% Dpt = Dc/DR*1e3; % pintle tip diameter (mm) - function of chamber diameter
Dpt = 25;
Ls = SR*Dpt; % skip length (mm) - distance outer flow travels before impingement point

% Pressure difference over injector
dP_ox = (Ps_ox-Pa)*1e5; % (Pa)
dP_f = (Ps_f-Pa)*1e5; % (Pa)

% Throttle
throttle = 1; % 1 = full throttle

%% Annular gap pintle

BF = 1; % Blockage factor, ratio of sum of circumference of orifices divided by circumference of pintle post
% Literature says 0.3 to 0.7 optimal, between 20 and 36 holes (ESA paper)

t_sleeve = 5.5; % thickness of sleeve (mm)
id_sleeve = Dpt - 2*t_sleeve; % sleeve ID (mm)

% Discharge coefficients for inner and outer flows, from experimental data https://www.researchgate.net/publication/301440576_Experiments_with_Pintle_Injector_Design_and_Development
Cd_i = 0.7; % MIT use 0.5 for cavitation https://wikis.mit.edu/confluence/display/RocketTeam/Modeling
Cd_o = 0.8; 

theta_pt = 40; % pintle tip angle (deg)
theta_post = theta_pt; % post angle (deg)
Dpr = 3; % pintle rod diameter (mm) ## Change this to be a dependent variable later
Dcg = 4.5; % center gap diameter (mm) ## Change this to be a dependent variable later
r_post = Dpt/2; % post diameter radius (mm)

A_o = mf_targ/((Cd_o*sqrt(2*rho_f*dP_f))*(1e-3)^2); % Outer orifice area (mm2); % Outer orifice area (mm2)
A_i = 25.5*throttle; % Inner orifice area (mm2) DETERMINED FROM DYER MODEL IN SEPARATE SCRIPT
Gap_o = sqrt(A_o/pi+(Dpt/2)^2)-Dpt/2; % Outer flow opening distance (mm)
Gap_i = A_i/(pi*id_sleeve*BF); % Pintle opening distance (mm)
Gap_iz = Gap_i/cosd(theta_pt); % Pintle axial opening distance (mm)

Dh_i = 2*Gap_i/1000; % Hydraulic diameter for inner flow (m)
Dh_o = 2*Gap_o/1000; % Hydraulic diameter for outer flow (m)

U_i = mox_targ/rho_ox/(A_i/(1e3)^2); % Velocity of inner flow (m/s)
U_o = mf_targ/rho_f/(A_o/(1e3)^2); % Velocity of outer flow (m/s)

mf_act_SPI = (Cd_o*sqrt(2*rho_f*dP_f))*A_o*(1e-3)^2; % Actual outer mass flow rate (kg/s), SPI model

Cd_eff_i = mox_targ/((A_i)*(1e-3)^2*sqrt(2*rho_f*Ps_f*1e5));
Cd_eff_o = mf_targ/((A_o)*(1e-3)^2*sqrt(2*rho_f*Ps_f*1e5));
Cd_eff = mp/((A_o+A_i)*(1e-3)^2*sqrt(2*rho_f*Ps_f*1e5));

%% Non-Dimensional Output Parameters

% Reynolds number
Re_i = rho_ox*U_i*Dh_i/mu_ox;
Re_o = rho_f*U_o*Dh_o/mu_f;

% Total momentum ratio
U_ia = U_i*sind((theta_pt+theta_post)/2); % Inner flow axial flow velocity
U_ir = U_i*cosd((theta_pt+theta_post)/2); % Inner flow radial flow velocity
mox_a = mox_targ*sind((theta_pt+theta_post)/2); % Inner flow axial mass flow
mox_r = mox_targ*cosd((theta_pt+theta_post)/2); % Inner flow radial mass flow
TMR = (mox_r*U_ir)/(mf_targ*U_o+mox_a*U_ia); % Total momentum ratio

% Momentum flux ratio
J = (rho_f*U_o^2)/(rho_ox*U_i^2);

% Inner and outer Weber numbers
We_i = rho_ox_g*U_i^2*Gap_i*1e-3/sigma_ox;
We_o = rho_f*U_o^2*Gap_o*1e-3/sigma_f;
% Different inner and outer We determines spray morphology
% Want We_i or We_o ~ 3000 for "fully developed fan spray"
% "The spray and atomization process and its effects on combustion performance of pintle injector,” Ph.D. thesis

% Weber number ### OXIDISER MUST BE GASEOUS AND FUEL LIQUID
% Uses relative velocity
U_rel = sqrt((U_ia-U_o)^2+U_ir^2);
We = rho_ox_g*U_rel^2*Gap_i*1e-3/sigma_f;

%% Spray Angle Prediction

% Half angle
theta_spray_half = acosd(1/(1+TMR));

%% Atomisation Prediction

zeta = (90-theta_pt)/(90);
q = 3.455-0.225*zeta;
SMD_Son = 1e3*Gap_i*zeta^(-1)*exp(4.0-q*(We)^0.1); % SMD in micron, Son https://doi.org/10.2514/6.2016-1453
SMD_Lee = -111.25*log(J^(-0.93)*We^0.05)+294.93; % SMD in micron, no consideration of pintle angle
