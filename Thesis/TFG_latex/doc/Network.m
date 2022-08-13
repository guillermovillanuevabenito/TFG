% DBS input onto a STN-GPe-GPi-Th network

% Rubin JE, Terman D (2004). High frequency stimulation of the subthalamic
% nucleus eliminates pathological thalamic rhythmicity in a computational
% model. J Comp Neurosci, 16:211-235

% Equations taken from ModelDB
% https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=116867&file=/rubin_terman2004/rubin_terman_pd.ode#tabs-1

% NORMAL STATE: 

%   Shift_th = -85 
%   Iapp_ge = 2
%   Gsyn_gege = 1 

% PSRKINSONIAN STATE

%   Shift_th = -80
%   Iapp_ge = -2.2
%   Gsyn_gege = 0

% These parameters must be changed in the function @f.

clearvars;
close all;
clc

tic

% Load initial conditions
load('ci.mat')
Tmax = 700;

% Solve equations
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,y] = ode45(@f,[0 Tmax],ini,options);

% Traces plot of TC
figure(1)
hold on
plot(t,y(:,1),'-b','linewidth',1);
axis([0 Tmax -110 0]);
set(gca,'fontsize',20);
xlabel('t  [ms]');
ylabel('V  [mV]');
title('Thalamic cell 1')

figure(2)
hold on
plot(t,y(:,2),'-b','linewidth',1);
axis([0 Tmax -110 0]);
set(gca,'fontsize',20);
xlabel('t  [ms]');
ylabel('V  [mV]');
title('Thalamic cell 2')

toc
