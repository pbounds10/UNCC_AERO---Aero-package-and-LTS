function [ vOut, accel ] = strLineVelo(dist, thrust, vIn, area, cd, mass )
%This calculates the straight line constant acceleration between points
%uses the thrust force and the aero forces to find acceleration
%Inputs:
%   dist - distance to travel (m)
%   thrust - thrust force (N)
%   vIn - initial velocity (m/s)
%   area - reference aero area (m^2)
%   cd - coefficient of drag ()
%   mass of vehicle (kg)
%Outputs:
%   vOut -  velocity final (m/s)
%   accel - accleration (m/s^2)

%Logic
%If the points are very close together then we should be able to
%approximate the velocity change between the points with simple kinmatics.
%To do this we are assuming that the distance is small enough that we can
%approximate the acceleration as constant.

den = 1.225; %kg/m^3 density of air

areoF = .5.*area.*den.*vIn.^2.*cd;

accel = (thrust-areoF)./mass;

vOut = sqrt( vIn.^2+2.*accel.*dist);



end

