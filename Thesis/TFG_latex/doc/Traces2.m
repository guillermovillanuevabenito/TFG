% Lambda-Omega Networks: the two-cell network
%
% This function computes x_{1}, y_{1}, x_{2} and y_{2} traces
%
% Output:
%   x1 : variable x_{1}
%   y1 : variable y_{1}
%   x2 : variable x_{2}
%   y2 : variable y_{2}

function [x1,y1,x2,y2] = Traces2(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha11,alpha12,alpha21,alpha22,dt,t)

x1 = zeros(1,length(t));
y1 = zeros(1,length(t));

% Initial conditions
theta1 = 0;  
x1(:,1) =sqrt(abs(lda1/b1))*cos(theta1);
y1(:,1) =sqrt(abs(lda1/b1))*sin(theta1);


x2 = zeros(1,length(t));
y2 = zeros(1,length(t));

x2(:,1) = sqrt(abs(lda2/b2))*cos(theta1);
y2(:,1) = sqrt(abs(lda2/b2))*sin(theta1);

% Connection
gmaxeff = 1;

for j=1:length(t)-1

    mod1 = x1(:,j).^2+y1(:,j).^2;

    k1x1 = lda1*x1(:,j) - omega1*y1(:,j) - (b1*x1(:,j)+a1*y1(:,j)).*mod1;
    k1x1 = k1x1 - c1*y1(:,j).*(1+x1(:,j)./sqrt(mod1))-(d1*y1(:,j).^3)./mod1 + gmaxeff*(alpha11.*x1(:,j)+ alpha12.*x2(:,j));

    k1y1 = omega1*x1(:,j) +lda1*y1(:,j) + (a1*x1(:,j)-b1*y1(:,j)).*mod1;    
    k1y1 = k1y1 + c1*x1(:,j).*(1+x1(:,j)./sqrt(mod1))+d1*(x1(:,j).*(y1(:,j).^2))./mod1;

    ax1 = x1(:,j)+k1x1*dt;
    ay1 = y1(:,j)+k1y1*dt;

    amod1 = ax1.^2+ay1.^2;


    mod2 = x2(:,j).^2+y2(:,j).^2;

    k2x1 = lda2*x2(:,j) - omega2*y2(:,j) - (b2*x2(:,j)+a2*y2(:,j)).*mod2;
    k2x1 = k2x1 - c2*y2(:,j).*(1+x2(:,j)./sqrt(mod2))-(d2*y2(:,j).^3)./mod2 + gmaxeff*(alpha22.*x2(:,j) + alpha21.*x1(:,j));

    k2y1 = omega2*x2(:,j) +lda2*y2(:,j) + (a2*x2(:,j)-b2*y2(:,j)).*mod2;    
    k2y1 = k2y1 + c2*x2(:,j).*(1+x2(:,j)./sqrt(mod2))+d2*(x2(:,j).*(y2(:,j).^2))./mod2;

    ax2 = x2(:,j)+k2x1*dt;
    ay2 = y2(:,j)+k2y1*dt;

    amod2 = ax2.^2+ay2.^2;


    k1x2 = lda1*ax1 - omega1*ay1 - (b1*ax1+a1*ay1).*amod1;
    k1x2 = k1x2 - c1*ay1.*(1+ax1./sqrt(amod1))-d1*ay1.^3./amod1 + gmaxeff*(alpha11.*ax1 + alpha12.*ax2);

    k1y2 = omega1*ax1 +lda1*ay1 + (a1*ax1-b1*ay1).*amod1;    
    k1y2 = k1y2 + c1*ax1.*(1+ax1./sqrt(amod1))+d1*(ax1.*(ay1.^2))./amod1; 
    
    x1(:,j+1) = x1(:,j)+(k1x1+k1x2)*dt/2;
    y1(:,j+1) = y1(:,j)+(k1y1+k1y2)*dt/2;

    k2x2 = lda2*ax2 - omega2*ay2 - (b2*ax2+a2*ay2).*amod2;
    k2x2 = k2x2 - c2*ay2.*(1+ax2./sqrt(amod2))-d2*ay2.^3./amod2 + gmaxeff*(alpha22.*ax2 + alpha21.*ax1);

    k2y2 = omega2*ax2 +lda2*ay2 + (a2*ax2-b2*ay2).*amod2;    
    k2y2 = k2y2 + c2*ax2.*(1+ax2./sqrt(amod2))+d2*(ax2.*(ay2.^2))./amod2;  

    x2(:,j+1) = x2(:,j)+(k2x1+k2x2)*dt/2;
    y2(:,j+1) = y2(:,j)+(k2y1+k2y2)*dt/2;
    
end

