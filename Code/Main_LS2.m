% Lambda-Omega Networks: the two-cell network
%
% This file compute 1-dimensional total-degenerated LSs on the connectivity
% parameter space. Symmetrical type-I heterogeneous network. Compensating parameter: alpha_{12}


clearvars;
close all;
clc
format long
 
tic

% Simulation parameters
Tmax = 100;
dt = 0.01;
t = 0:dt:Tmax;

tmin = 0.75*Tmax;
tmax = Tmax;

% Compensating parameter
alpha12 = 1:0.1:3;

% Parameters cell 1
a1 = 1;
lda1 = 1;
b1 = 1;

c1 = 0;
d1 = 0;

omega1 = 1;

% Parameters for cell 2
a2 = 1;
lda2 = 3;
b2 = 3;

c2 = c1;
d2 = d1;

omega2 = 1;

% Vectors to store solution
M1 = zeros(length(alpha12),1);
M2 = zeros(length(alpha12),1);
M3 = zeros(length(alpha12),1);

% LS Attributes
f0 = 0.3868;
amp0 = 1.5;

% LS parameters
eps = 0.5;
niter = 15;

% Initial value
alpha21_aprox = 1;

for i1 = 1:length(alpha12)

    % Target function
    f = @(x) f_target(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha12(i1),x,dt,t,tmin,tmax,amp0)-f0;

    % Bisection-like algorithm
    [alpha21_aprox, alpha11_aprox, alpha22_aprox, error]= FindLevelSet(alpha21_aprox,eps,f,niter,lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha12(i1),dt,t,tmin,tmax,amp0,f0);

    %Update
    M1(i1) = alpha21_aprox;
    M2(i1) = alpha11_aprox;
    M3(i1) = alpha22_aprox;

    % Check error
    if error > 0.01
        break
    end

end

% Plots
figure(1)
hold on
plot(alpha12,M1,'-b','linewidth',2);
set(gca,'fontsize',20);
xlabel('\alpha_{12}');
ylabel('\alpha_{21}');

figure(2)
plot(alpha12,M2,'-r','linewidth',2);
set(gca,'fontsize',20);
xlabel('\alpha_{12}');
ylabel('\alpha_{11}');

figure(3)
hold on
plot(alpha12,M3,'-g','linewidth',2);
set(gca,'fontsize',20);
xlabel('\alpha_{12}');
ylabel('\alpha_{22}');


toc