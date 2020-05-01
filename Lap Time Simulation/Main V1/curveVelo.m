function [ vOut ] = curveVelo( mass, rad, area, cd, tireCo, cdf, vIn )
%calculates the max velocity in the curve
%   Detailed explanation goes here

den = 1.225; %kg/m^3
g = 9.81;

normForce = mass*9.81 + cdf*.5*vIn.^2*den*area;

a = normForce*tireCo;

% b = mass.^2./(rad.^2);
% 
% c = (.5*den*area*cd)^2;
% 
% d = a./(b+c);
% 
% vOut = d^(.25);

vOut = ( (a*rad)/mass )^.5;


end

