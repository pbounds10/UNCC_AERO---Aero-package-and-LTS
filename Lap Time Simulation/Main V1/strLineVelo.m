function [ vOut ] = strLineVelo( dist, powEff, vIn, area, cd, mass )
%This calculates the straight line constant acceleration
%   Detailed explanation goes here

den = 1.225; %kg/m^3 density of air

num = (powEff./vIn)-(.5.*vIn.^2*area*cd*den);
num = num .* dist .*2;

den = mass;

v2 = num./den + vIn.^2;

vOut = v2.^(.5);



end

