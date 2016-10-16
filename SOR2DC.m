function [ u,v ] = SOR2DC( sigma, x, y, xj, yj, x2j, y2j )
if (x == (xj + x2j)/2 && y == (yj + y2j)/2 )
    u = 0;
    v=sigma * 0.5;
else
    u = (sigma / (4 * pi)) * log( ((x - xj)^2 + (y - yj)^2) / ((x - x2j)^2 + (y - y2j)^2) );
    v = (sigma / (2 * pi)) * (atan((y - y2j) / (x - x2j)) - atan((y - yj) / (x - xj)) );
end
end

