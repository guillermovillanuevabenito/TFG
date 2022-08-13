% Auxiliary function
%
% For a given pair of cros-connectivity parameter this function computes
% self-connectivity parameters alpha_{11} and alpha_{22} preserving a
% desired network amplitue (amp0)
%
% Output
%   alpha11aprox : value of compensated parameter alpha_{11}
%   alpha22aprox : value of compensated parameter alpha_{22}

function [alpha11aprox,alpha22aprox] = Selfpar(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha12,alpha21,dt,t,tmin,tmax,amp0)

q = 2;

alpha11aprox = 0;
alpha22aprox = 0;
    
alpha22 = -alpha21-3:1:-alpha21+6.5; 
alpha11 = -alpha12-3:1:-alpha12+6.5; 

p = 1;
pp = 9;

% Search-method to compute compensated parameters alpha_{11} and \alpha_{22}
for i1 = 1:pp

    if i1>1

        alpha11left = alpha11aprox-p/(q^(i1-2));
        alpha11right = alpha11aprox + p/(q^(i1-2));
        alpha22left = alpha22aprox-p/(q^(i1-2));
        alpha22right = alpha22aprox + p/(q^(i1-2));
        alpha11 = alpha11left:p/q^(i1-1):alpha11right;
        alpha22 = alpha22left:p/q^(i1-1):alpha22right;

    end

    [alpha11aprox,alpha22aprox] = Iter(lda1,b1,omega1,a1,c1,d1,lda2,b2,omega2,a2,c2,d2,alpha11,alpha12,alpha21,alpha22,dt,t,tmin,tmax,amp0);
    alpha11aprox = alpha11aprox(end);
    alpha22aprox = alpha22aprox(end);

end

