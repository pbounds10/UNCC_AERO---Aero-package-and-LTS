function [ vOut ] = curveAccelVelo( mass, rad, area, cd, tireCo, cdf, vIn, dist  )
%This calculates the velocity of the car in a curve as it accelerates
%   Detailed explanation goes here

den = 1.225; %kg/m^3
g = 9.81;

ac = vIn^2/rad;


normForce = mass*9.81 + cdf*.5*vIn.^2*den*area;

tireForce = normForce*tireCo;

vMax = ( (tireForce*rad)/mass )^.5;

vOut = (vIn^2+2*ac*dist)^(.5);

if (vOut > vMax)
    vOut = vMax;

end


end

