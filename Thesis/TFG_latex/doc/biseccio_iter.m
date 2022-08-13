% Bisection method

function [vect_x,vect_r,pos] = biseccio_iter(Ini,niter,f)

    vect_x = zeros(1,niter); 
    vect_r = zeros(1,niter); 
    x0 = Ini(1); f0 = f(x0); 
    a =  Ini(2); fa = f(a); 
    
    if f0*fa > 0 
        error('Interval inicial inadequat') 
    end
    
    ok = 0;
    i = 1;
    pos = niter;
    tol = 1e-5;
    
    while (i <= niter && ok == 0)
        x1 = (x0 + a)/2;
        f1 = f(x1);
        vect_x(i) = x0;  
        vect_r(i) = abs((x1-x0)/x1);
        
        if abs(f(x0))< tol
            pos = i;
            ok = 1;
        else
            if (f1*f0 < 0)
                a = x0; 
            end
            x0 = x1; f0 = f1; 
            i = i+1;
        end
    end