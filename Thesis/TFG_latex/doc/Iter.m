% Auxiliary function
%
% Find the minimum of the function sqrt((amp1-amp0)^2 + (amp2-amp0)^2) on
% an given region on the alpha11-alpha22 parameter space where amp1 is the
% amplitude of cell-1, amp2 is the amplitude os cell-2 and amp0 is the
% amplitude value defining the desired LS on connectivity parameter space
%
% Output:
%   sol1: value of parameter alpha_{11}
%   sol2: value of parameter alpha_{22}

function [sol1,sol2] = Iter(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha11,alpha12,alpha21,alpha22,dt,t,tmin,tmax,amp0)

ampn11 = zeros(length(alpha11),length(alpha22));
ampn22 = zeros(length(alpha11),length(alpha22));
ampdif = zeros(length(alpha11),length(alpha22));


for i1 = 1:length(alpha11)

    for i2 = 1:length(alpha22)

        [x1,~,x2,~] = Traces2(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha11(i1),alpha12,alpha21,alpha22(i2),dt,t);

        [amp1,~,osc1] = Oscillation(x1,tmin,tmax,t,dt);
        [amp2,~,osc2] = Oscillation(x2,tmin,tmax,t,dt);

        if osc1 == 1 && osc2 == 1

             ampdif(i1,i2) = sqrt((amp1-amp0)^2 + (amp2-amp0)^2);
             
             ampn22(i1,i2) = amp2;
             ampn11(i1,i2) = amp1;

        else
             ampdif(i1,i2) = 10;
             ampn22(i1,i2) = 10;
             ampn11(i1,i2) = 10;
        end

    end

end
[x,y] = find(ampdif==min(ampdif(:)));
sol1 = alpha11(x);
sol2 = alpha22(y);