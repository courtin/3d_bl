function [mids,dx] = midpoints(x)
    %Places points at the midpoint of an array of x,y,z points.
    S = size(x);
    N = S(2);
    mids = zeros(S(1), N-1);
    dx = zeros(S(1), N-1);
    for i = 1:N-1
        mids(:,i) = (x(:,i+1)+x(:,i))./2;
        dx(:,i) = (x(:,i+1)-x(:,i));
    end

end