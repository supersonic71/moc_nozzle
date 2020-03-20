%this function expects input in degrees
function [xo,yo] = slope_inter(X, Y, m1, m2)
%no way to take [x1, x2], [y1, y2] as input?
    m1 = tand(m1);
    m2 = tand(m2);


    x1 = X(1); x2 = X(2); y1 = Y(1); y2 = Y(2);
    c1 = y1-m1*x1; c2 = y2-m2*x2;
    xo = (c2-c1)./(m1-m2);
    yo = m1*xo + c1;
%line([x1 xo],[y1 yo]);
%line([x2 xo],[y2 yo]);
%axis([0 6 0 6]);
end