function [V] = vorvel2(r, ra, rb, LBOUND)
%Calculate velocity at a point r due to a bound vortex defined by trailing
%legs starting at ra and rb. 

b = ra-r;
a = rb-r;

xh = [1;0;0];

na = norm(a);
nb = norm(b);

if LBOUND == 0
    T1 = 0;
else
    T1 = cross(a,b)/(na*nb+dot(a,b))*(1/na+1/nb);
end
T2 = cross(a,xh)/(na-dot(a,xh))*(1/na);
T3 = cross(b,xh)/(nb-dot(b,xh))*(1/nb);

V = 1/(4*pi)*(T1 + T2 - T3);
end

