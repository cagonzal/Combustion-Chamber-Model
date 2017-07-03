function [ K ] = K_methane( T, P, phi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

gas = GRI30('Multi');
fuel = GRI30('Multi');
ox = GRI30('Multi');
FOst = 0.25; % Stoichiometric Fuel-Oxidizer Ratio
FO = phi * FOst;
R = 8.314;
nu = 1/phi/FOst; % OF Ratio
m = 16; % Molecular weight of methane
rho = 300; % Value Pablo says to use for liquid methane
m_pr = 80; % Molecular weight of products of reaction

nsp = nSpecies(gas);
iLOX = speciesIndex(gas,'O2');
iCH4 = speciesIndex(gas,'CH4');

y = zeros(nsp,1);
y(iCH4,1) = FO / (FO+1);
y(iLOX,1) = 1 / (FO + 1);

% Atmospheric Boiling Conditions for Methane
t_b0 = 109;
P0 = 101325;

% Assume these are constant, taken at STP
hfg = 509;
hc = 50016;

% Use Clausius-Clapeyron Equation to get boiling temp at our pressure
A = P0;
B = R / hfg / m;
tb = (1/t_b0 - B*log(P/A))^(-1);

% Set state for gas using tbar and p
set(gas,'Temperature',T,'Pressure',P,'Y',y);

equilibrate(gas,'HP',1);
tad = temperature(gas);

% Guess that Tf = adiabatic flame temperature
tf = tad;
tbar = 0.5*(tb + tf);
%Tinf = 300;

set(fuel,'Temperature',tbar,'Pressure',P,'MoleFractions','CH4:1');
set(ox,'Temperature',tbar,'Pressure',P,'MoleFractions','O2:1');

kf = thermalConductivity(fuel);
kox = thermalConductivity(ox);
k = 0.4*kf + 0.6*kox; % Emprical relation for thermal conductivity
cp = cp_mass(fuel);

Boq = (hc/nu * 1000 + cp*(T-tb))/(hfg * 1000);
K = 8*k/(rho*cp)*log(1+Boq);

end
