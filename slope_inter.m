function [xo,yo] = slope_inter(X, Y, m1, m2)

% The function returns the intersecting point of two lines.
% A line is defined by a point and a corresponding slope
% Two points and two respective slopes are inputs
% x coordinates of points go in X, y coordinates of points go in Y
% The slopes are defined in degrees 

    m1 = tand(m1); %Converting degrees to actual slope values
    m2 = tand(m2);

    % TODO : no way to take [x1, x2], [y1, y2] as input?
    
    % Decomposing X and Y into respective points
    x1 = X(1); x2 = X(2); y1 = Y(1); y2 = Y(2);
    
    % Formulas from basic line geometry
    
    c1 = y1-m1*x1;
    c2 = y2-m2*x2;
    
    xo = (c2-c1)./(m1-m2);
    yo = m1*xo + c1;
    
    % Test visualisation
    % line([x1 xo],[y1 yo]);
    % line([x2 xo],[y2 yo]);
    
end