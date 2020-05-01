function [ fFront,fRear] = brakeByPress( pressIn, presMax, wfs, wrs, radWheelF,radWheelR, bgf, bgr, pCrit, pro)
%UNTITLED Summary of this function goes here
%INPUTS:
%   pressIn - pressure into system in kpa
%   presMax - maximum braking pressure in kpa
%   wfs/wrs - static front and rear weight in N
%   radWheel - radius of wheel in m
%   bgf/bgr - brake gain front and rear in N-m/kpa
%   pCrit - critical brake pressure in kpa
%   pro - proportioning
%   g - gravitational Acceleration in m/s^2
%   cgH - center of gravity height in m
%   len - wheelbase in m
%   lenA - distance from cg to front axle - m

g = 9.81;

weight = wfs + wrs; %in N
mass = weight./g; %in kg

if ( pressIn > presMax)
    pressIn = presMax;
end

%Check to see if the brake proportioning has engaged
if (pressIn > pro)
    pf = pressIn;
    pr = pressIn;
else
    pf = pressIn;
    pr = pCrit + pro.*(pressIn-pCrit);
end

fFront = bgf.*pf./radWheelF; %front braking force
fRear = bgr.*pr./radWheelR; %rear braking force

dx = ( fFront + fRear )./weight; %acceleration in g's
dx = dx.*g; %in m/s^2

brakeForce = fFront+fRear; %N


end

