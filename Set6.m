%% Set Constants and Initial Conditions
clc;
clear all; 
close all;

global P0 A mdot_ox mdot_l0 rho_l FO_st hfg ii

R = 8.314;
m_ox = 32;
m_f = 16;
vd0 = 10;
T0 = 600;
P0 = 3.4474 * 10^(6);
D0 = 10^(-6) * [30 50 80 100 200];
A = 0.157;
Ainj = 0.0157;
L = 0.75;
phi0 = 1.139;
phi_g0 = 0.45;
FO_st = m_f / (2*m_ox);
FO_0 = phi0 * FO_st;
FO_0g = phi_g0 * FO_st;
hfg = 509*1000;

rho_l = 300; % In Turns, given to us in email from Pablo
mdot_l0 = rho_l*Ainj*vd0;
mdot_fg0 = mdot_l0*phi_g0 / (phi0-phi_g0);
mdot_ox = mdot_fg0 * 4 / phi_g0;
mdot_g0 = mdot_fg0 + mdot_ox;
%mdot_total = mdot_l0 + mdot_g0;

gas = GRI30;
nsp = nSpecies(gas);
iLOX = speciesIndex(gas,'O2');
iCH4 = speciesIndex(gas,'CH4');

y_0 = zeros(nsp,1);
y_0(iCH4,1) = FO_0 / (1 + FO_0);
y_0(iLOX,1) = 1 / (1 + FO_0);

set(gas,'Temperature',T0,'Pressure',P0,'Y',y_0);
equilibrate(gas,'HP');
Tg0 = temperature(gas);

%% Use ODE45 to calculate everything

% Make sure to save all the data

%xspan = linspace(0,0.75);
%ii = 1;
%[x1,z1] = ode45(@ode_rhs, xspan, [D0(1)^2, Tg0, vd0, mdot_g0, phi_g0]);
%save('set6_d1','x1','z1');
%ii = 1;
%[x2,z2] = ode45(@ode_rhs, xspan, [D0(2)^2, Tg0, vd0, mdot_g0, phi_g0]);
%ii = 1;
%save('set6_d2','x2','z2');
%[x3,z3] = ode45(@ode_rhs, xspan, [D0(3)^2, Tg0, vd0, mdot_g0, phi_g0]);
%ii = 1;
%save('set6_d3','x3','z3');
%[x4,z4] = ode45(@ode_rhs, xspan, [D0(4)^2, Tg0, vd0, mdot_g0, phi_g0]);
%ii = 1;
%save('set6_d4','x4','z4');
%[x5,z5] = ode45(@ode_rhs, xspan, [D0(5)^2, Tg0, vd0, mdot_g0, phi_g0]);
%save('set6_d5','x5','z5');

%% Load Data

%ode45 takes very long to run so we load in previously saved data

load('set6_d1.mat');
load('set6_d2.mat');
load('set6_d3.mat');
load('set6_d4.mat');
load('set6_d5.mat');

%% Calculate stuff not given out by oderhs

for jj = 1:length(x1)
    gas1 = GRI30('Multi');
    
    Tg1 = z1(jj,2);
    Tg2 = z2(jj,2);
    Tg3 = z3(jj,2);
    Tg4 = z4(jj,2);
    Tg5 = z5(jj,2);
    
    phi_g1 = z1(jj,5);
    phi_g2 = z2(jj,5);
    phi_g3 = z3(jj,5);
    phi_g4 = z4(jj,5);
    phi_g5 = z5(jj,5);
    
    FO1 = phi_g1 * FO_st;
    FO2 = phi_g1 * FO_st;
    FO3 = phi_g1 * FO_st;
    FO4 = phi_g1 * FO_st;
    FO5 = phi_g1 * FO_st;
    
    y1 = zeros(nsp,1);
    y1(iCH4,1) = FO1 / (1 + FO1);
    y1(iLOX,1) = 1 / (1 + FO1);
  
    y2 = zeros(nsp,1);
    y2(iCH4,1) = FO2 / (1 + FO2);
    y2(iLOX,1) = 1 / (1 + FO2);
    
    y3 = zeros(nsp,1);
    y3(iCH4,1) = FO3 / (1 + FO3);
    y3(iLOX,1) = 1 / (1 + FO3);
    
    y4 = zeros(nsp,1);
    y4(iCH4,1) = FO4 / (1 + FO4);
    y4(iLOX,1) = 1 / (1 + FO4);
    
    y5 = zeros(nsp,1);
    y5(iCH4,1) = FO5 / (1 + FO5);
    y5(iLOX,1) = 1 / (1 + FO5);
    
    set(gas1,'Temperature',Tg1,'Pressure',P0,'Y',y1);
    rho_g1 = density(gas1);
    set(gas1,'Temperature',Tg2,'Pressure',P0,'Y',y2);
    rho_g2 = density(gas1);
    set(gas1,'Temperature',Tg3,'Pressure',P0,'Y',y3);
    rho_g3 = density(gas1);
    set(gas1,'Temperature',Tg4,'Pressure',P0,'Y',y4);
    rho_g4 = density(gas1);
    set(gas1,'Temperature',Tg5,'Pressure',P0,'Y',y5);
    rho_g5 = density(gas1);
end

vg1 = z1(:,4) / A / rho_g1;
vg2 = z2(:,4) / A / rho_g2;
vg3 = z3(:,4) / A / rho_g3;
vg4 = z4(:,4) / A / rho_g4;
vg5 = z5(:,4) / A / rho_g5;

mdot_l1 = z1(:,3) / Ainj / rho_l;
mdot_l2 = z2(:,3) / Ainj / rho_l;
mdot_l3 = z3(:,3) / Ainj / rho_l;
mdot_l4 = z4(:,3) / Ainj / rho_l;
mdot_l5 = z5(:,3) / Ainj / rho_l;

%% Generate Plots

figure(1);
set(gcf,'units','normalized','position',[0 0 1 1])
hold on;
plot(x1,sqrt(z1(:,1))/D0(1));
plot(x2,sqrt(z2(:,1))/D0(2));
plot(x3,sqrt(z3(:,1))/D0(3));
plot(x4,sqrt(z4(:,1))/D0(4));
plot(x5,sqrt(z5(:,1))/D0(5));
title('Normalized Droplet Diameter');
xlabel('Distance from Inlet [m]');
ylabel('D/D0')
legend('D0 = 30 \mum','D0 = 50 \mum','D0 = 80 \mum','D0 = 100 \mum',...
    'D0 = 200 \mum')

figure(2);
set(gcf,'units','normalized','position',[0 0 1 1])
hold on;
plot(x1,z1(:,2));
plot(x2,z2(:,2));
plot(x3,z3(:,2));
plot(x4,z4(:,2));
plot(x5,z5(:,2));
title('Droplet Temperature');
xlabel('Distance from Inlet [m]');
ylabel('Tg [K]')
legend('D0 = 30 \mum','D0 = 50 \mum','D0 = 80 \mum','D0 = 100 \mum',...
    'D0 = 200 \mum')


figure(3);
set(gcf,'units','normalized','position',[0 0 1 1])
hold on;
plot(x1,z1(:,3));
plot(x2,z2(:,3));
plot(x3,z3(:,3));
plot(x4,z4(:,3));
plot(x5,z5(:,3));
title('Droplet Velocity');
xlabel('Distance from Inlet [m]');
ylabel('Vg [m/s]')
legend('D0 = 30 \mum','D0 = 50 \mum','D0 = 80 \mum','D0 = 100 \mum',...
    'D0 = 200 \mum')

figure(4);
set(gcf,'units','normalized','position',[0 0 1 1])
hold on;
plot(x1,vg1);
plot(x2,vg2);
plot(x3,vg3);
plot(x4,vg4);
plot(x5,vg5);
title('Gas Velocity');
xlabel('Distance from Inlet [m]');
ylabel('Vg [m/s]')
legend('D0 = 30 \mum','D0 = 50 \mum','D0 = 80 \mum','D0 = 100 \mum',...
    'D0 = 200 \mum')

figure(5);
set(gcf,'units','normalized','position',[0 0 1 1])
hold on;
plot(x1,z1(:,4));
plot(x2,z2(:,4));
plot(x3,z3(:,4));
plot(x4,z4(:,4));
plot(x5,z5(:,4));
title('Gas Mass Flow Rate');
xlabel('Distance from Inlet [m]');
ylabel('Mdot_g [kg/s]')
legend('D0 = 30 \mum','D0 = 50 \mum','D0 = 80 \mum','D0 = 100 \mum',...
    'D0 = 200 \mum')

figure(6);
set(gcf,'units','normalized','position',[0 0 1 1])
hold on;
plot(x1,mdot_l1);
plot(x2,mdot_l2);
plot(x3,mdot_l3);
plot(x4,mdot_l4);
plot(x5,mdot_l5);
title('Liquid Mass Flow Rate');
xlabel('Distance from Inlet [m]');
ylabel('Mdot_g [kg/s]')
legend('D0 = 30 \mum','D0 = 50 \mum','D0 = 80 \mum','D0 = 100 \mum',...
    'D0 = 200 \mum')

figure(7);
set(gcf,'units','normalized','position',[0 0 1 1])
hold on;
plot(x1,z1(:,5));
plot(x2,z2(:,5));
plot(x3,z3(:,5));
plot(x4,z4(:,5));
plot(x5,z5(:,5));
title('Equivalence Ratio');
xlabel('Distance from Inlet [m]');
ylabel('phi_g')
legend('D0 = 30 \mum','D0 = 50 \mum','D0 = 80 \mum','D0 = 100 \mum',...
    'D0 = 200 \mum')

%% Effect of Injection Velocity

%xspan = linspace(0,0.75);

%vd0 = 20;
%mdot_l0 = rho_l*Ainj*vd0;
%ii = 1;
%[x20,z20] = ode45(@ode_rhs, xspan, [D0(2)^2, Tg0, vd0, mdot_g0, phi_g0]);
%save('vd0s_20','x20','z20');

%vd0 = 30;
%mdot_l0 = rho_l*Ainj*vd0;
%ii = 1;
%[x30,z30] = ode45(@ode_rhs, xspan, [D0(2)^2, Tg0, vd0, mdot_g0, phi_g0]);
%save('vd0s_30','x30','z30');

%vd0 = 50;
%mdot_l0 = rho_l*Ainj*vd0;
%ii = 1;
%[x50,z50] = ode45(@ode_rhs, xspan, [D0(2)^2, Tg0, vd0, mdot_g0, phi_g0]);
%save('vd0s_50','x50','z50');

%% Load Injection Velocity Data

load('vd0s_20.mat');
load('vd0s_30.mat');


%% Plot effect of injection velocity

figure(8);
set(gcf,'units','normalized','position',[0 0 1 1])
hold on;
plot(x1,sqrt(z2(:,1))/D0(2));
plot(x20,sqrt(z20(:,1))/D0(2));
plot(x30,sqrt(z30(:,1))/D0(2));
title('Normalized Droplet Diameter (D0 = 20 \mum)');
xlabel('Distance from inlet [m]');
ylabel('D/D0');
legend('v0 = 10 m/s','v0 = 20 m/s','v0 = 30 m/s');
