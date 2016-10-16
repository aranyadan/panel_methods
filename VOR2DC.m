function [ u,v ] = VOR2DC( gamma, x, y, xj, yj, x2j, y2j )
%% function to calculate velocity due to panel vortexes
if (x == (xj + x2j)/2 && y == (yj + y2j)/2 )
    u = gamma * -0.5;
    v=0;
else
    u = (gamma / (2 * pi)) * (atan((y - y2j) / (x - x2j)) - atan((y - yj) / (x - xj)) );
    v = (-1 * gamma / (4 * pi)) * log( ((x - xj)^2 + (y - yj)^2) / ((x - x2j)^2 + (y - y2j)^2) );
end
end

