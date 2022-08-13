% Lambda-Omega Networks: the self-connected cell
% 
% This function compute x and y traces
%
% Output:
%   x = variable x
%   y = variable y

function [x,y] = Traces(lda,b,omega,a,c,d,alpha,dt,t)

% Parameters
n = length(lda);
x = zeros(length(t),n);
y = zeros(length(t),n);

% Initial conditions
x(1,:) = 1;
y(1,:) = 0;

% Self-connection
gmaxeff = 1;

% Simulation
for j=1:length(t)-1

    mod = x(j,:).^2+y(j,:).^2;

    kx1 = lda.*x(j,:) - omega.*y(j,:) - (b.*x(j,:)+a*y(j,:)).*mod;
    kx1 = kx1 - c*y(j,:).*(1+x(j,:)./sqrt(mod))-(d*y(j,:).^3)./mod + gmaxeff*alpha.*x(j,:);

    ky1 = omega.*x(j,:) +lda.*y(j,:) + (a*x(j,:)-b.*y(j,:)).*mod;    
    ky1 = ky1 + c*x(j,:).*(1+x(j,:)./sqrt(mod))+d*(x(j,:).*(y(j,:).^2))./mod;

    ax = x(j,:)+kx1*dt;
    ay = y(j,:)+ky1*dt;

    amod = ax.^2+ay.^2;

    kx2 = lda.*ax - omega.*ay - (b.*ax+a*ay).*amod;
    kx2 = kx2 - c*ay.*(1+ax./sqrt(amod))-d*ay.^3./amod + gmaxeff*alpha.*ax;

    ky2 = omega.*ax +lda.*ay + (a*ax-b.*ay).*amod;    
    ky2 = ky2 + c*ax.*(1+ax./sqrt(amod))+d*(ax.*(ay.^2))./amod;  

    x(j+1,:) = x(j,:)+(kx1+kx2)*dt/2;
    y(j+1,:) = y(j,:)+(ky1+ky2)*dt/2;
    
end