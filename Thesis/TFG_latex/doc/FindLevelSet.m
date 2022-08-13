% Lambda-Omega networks: the two-cell network
%
% Auxiliary funtion: compute compensated connectivity parameters for a
% given value of parameter alpha_{12} (compensating parameter). For each
% parameter alpha_{21}, self-connectivity parameter alpha_{11} and
% alpha_{22} are computed to preserve the desired amplitude (amp0) in order
% to a have a 1-d problem which consist on finding parameter alpha21
% preserving the desired frequency (f0).
%
% Output:
%   alpha21_aprox : compensated parameter alpha{21}
%   alpha11_aprox : compensated parameter alpha_{11}
%   alpha22_aprox : compensated parameter alpha_{22}
%   error : sqrt((amp1-amp0)^2 + (amp2-amp0)^2 + (f-f0)^2 where amp1 is the
%   amplitude of cell-1, amp2 is the amplitude of cell-2, f is the network
%   frequency, and amp0 and f0 define the desired LS on connectivity
%   parameter space


function [alpha21_aprox, alpha11_aprox, alpha22_aprox, error] = FindLevelSet(ini,eps,f,niter,lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha12,dt,t,tmin,tmax,amp0,f0)


f1 = f(ini+eps);
f2 = f(ini-eps);
ff0 = f(ini);
ok = 0;

% Bisection-like method to the target function f
if (sign(f1*ff0) == -1)

    [vect_x,~,p] = biseccio_iter([ini ini+eps],niter,f);
    alpha21_aprox = vect_x(p);

    ok = 1;

elseif (sign(f2*ff0) == -1)

    [vect_x,~,p] = biseccio_iter([ini-eps ini],niter,f);
    alpha21_aprox = vect_x(p);

    ok = 1;

end

if (ok == 1)
    % From the aproximation of compensated parameter we compute de error
    [alpha11_aprox,alpha22_aprox] = Selfpar(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha12,alpha21_aprox,dt,t,tmin,tmax,amp0);

    [x1,~,x2,~] = Traces2(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha11_aprox,alpha12,alpha21_aprox,alpha22_aprox,dt,t);

    [amp1,f,~] = Oscillation(x1,tmin,tmax,t,dt);
    [amp2,~,~] = Oscillation(x2,tmin,tmax,t,dt);

    error = sqrt((amp1-amp0)^2 + (amp2-amp0)^2 + (f-f0)^2);
    
    X = ['ERROR: ', num2str(error),';  alpha21: ', num2str(alpha21_aprox),';  alpha11: ', num2str(alpha11_aprox),';  alpha22: ', num2str(alpha22_aprox) ];
    disp(X)
   
else
    
    X = ('¡¡¡ WARNING !!!: Zero not found. The value of eps might be changed');
    disp(X)
   
end

