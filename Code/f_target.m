% Target function
%
% Output
%   f : frequency on the the point of the amplitude LS on connectivity
%   parameters space (A_{Net} = amp0) defined by the cross-connectivity
%   parameters alpha_{12} and \alpha_{21}

function [f,alpha11aprox,alpha22aprox] = f_target(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha12,alpha21,dt,t,tmin,tmax,amp0)

[alpha11aprox,alpha22aprox] = Selfpar(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha12,alpha21,dt,t,tmin,tmax,amp0);


[x1,~,~,~] = Traces2(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha11aprox,alpha12,alpha21,alpha22aprox,dt,t);

[~,f,~] = Oscillation(x1,tmin,tmax,t,dt);

