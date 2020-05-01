function [ vOut ] = constABraking( accel, vIn, dist )
%brakes the car using a constant negtive acceleration
%   Detailed explanation goes here

vOut = ( vIn.^2+2*accel*dist )^.5; %This should increase speed as you are trying to where the straight speed at a point can decel to make the corner speed

end

