% This function compute the amplitude and frequency of a given
% sinusoidal-like voltage signal
%
% Output:
%   amp = Amplitude
%   f = Frequency
%   osc = 1 if the amplitude and frequency are well computed or 0 otherwise

function [amp,f,osc] = Oscillation(V,tmin,tmax,t,dt)

jmin = floor(tmin/dt);
jmax = floor(tmax/dt);

Vmax = zeros(1,1);
Vmin = zeros(1,1);
peakt = zeros(1,1);

cnt1 = 0;
cnt2 = 0;

for j=jmin+2:jmax-1 
    
    if V(j)>V(j-1) && V(j)>V(j+1)
        cnt1=cnt1+1;
        Vmax(cnt1)=V(j);
        peakt(cnt1)=t(j);
    end
    
    if V(j)<V(j-1) && V(j)<V(j+1)
        cnt2=cnt2+1;
        Vmin(cnt2) = V(j);
    end
    
end

Vmax = flip(Vmax);
Vmin = flip(Vmin);
peakt = diff(peakt);

T = flip(peakt);
e = diff(Vmax);
e = abs(e);

% Avoid non-sustained oscillations
if isempty(e)
    
    amp = 0;
    f = 0;
    osc = 0;
    
else

    aux = (Vmax(1)-Vmin(1))/2;
  
    % Avoid damped and non-sinusoidal-like oscillations
    if (e(1) < 0.01 && aux > 0.1)
        amp = (mean(Vmax)-mean(Vmin))/2;
        f = 1/mean(T);
        osc = 1;
    else
        amp = 0;
        f = 0;
        osc = 0;
    end
end
