function [ u,v ] = VOR2D( gamma, x, y, xj, yj )
%% function to calculate velocity due to discrete vortexes
r = sqrt( (x - xj)^2 + (y - yj)^2);
pos = [(x-xj) ; (y-yj)];
unit = [0 1 ; -1 0];
velocity = (gamma / (2 * pi * r * r)) * unit * pos;
u=velocity(1);
v=velocity(2);
end

