function [ax,wfs,wrs,kappaPeak] = longAccel(u,v, gRatio, fd, torque, wheelRad, rpmArray, wfs, wrs, cgH, len,lenA, SA, tireCos, camber, driveWheels, springF,springR, yaw, pitch, aeroCD, aeroCDF, aeroFront, area, density, symetric, curveRad)
%UNTITLED Summary of this function goes here
%Inputs:
%   wfs - static weight of front
%   wrs - static weight of rear
%   cgH - height of the center of gravity in m
%   len - length of the car in m
%   accel - acceleration of car m/s^2
%   lenA - distance from front axle to cg
%   u - long velocity m/s
%   v = lat velocity m/s
%   curveRad - radius of curve in m if straight use NaN
%   gRatio - array of gearRatios
%   fd - final drive
%   torque - array of torques that are from dyno data
%   wheelRad - radius of drive wheels
%   rpmArray - rpms at which dyno torque occured
%   driveWheels - rwd = 1,fwd = 2
%   SA - slip angle
%   camber - camber angle in degrees of drive wheels
%   co - the load influenced and camber influece tire coefficents
%   output - 1 = laterial, 3 = corrective moment, and 2 = longitudinal,
%   kappa - slip ratio 0-1
%Output:

%% Initialize
mass = (wfs+wrs)./9.81; %kg
velo = sqrt(u.^2+v^2);
wfs_neutral = wfs;
wrs_neutral = wrs;

%% 1) find thrust made my engine
    thrust = thrustForce1(velo, gRatio, fd, torque, wheelRad, rpmArray); %thrust force in N produced by engine at wheels

%% Loop til steady
while(1) 
    %% 2) find max traction
    if (driveWheels == 1) %rwd
        [~,~,~,kappaPeak] = magicTire(wrs./2000,0,0,camber,tireCos,2);
        tireForce = 2.*magicTire(wrs./2000,SA,kappaPeak,camber,tireCos,2);
    else
        [~,~,~,kappaPeak] = magicTire(wfs/2000,0,0,camber,tireCos,2);
        tireForce = 2.*magicTire(wfs/2000,SA,kappaPeak,camber,tireCos,2);
    end

    %% 3) find limit of force tires or engine
    if (thrust < tireForce) %enine limited
        Fx = thrust;
    else %tire limited
        Fx = tireForce;
        %disp('tires')
    end

    %% 4) find acceleration
    drag = aeroForce(yaw,pitch,u,v,aeroCD,aeroCDF,aeroFront,area,density,1);
    if (~isnan(curveRad))
        ax = (Fx-drag)./mass -(v.*u)./curveRad;
    else
        ax = (Fx-drag)./mass;
    end

    %% 5) find weight transfer
    wfsOld = wfs;
    [wfs,wrs] = weightTrans(wfs, wrs, cgH, len,lenA, ax );
    
    %% 6) find ride height change and pitch
        %currently not do but could be

    %% 7) check to see if its time to break loop
    if (abs(wfsOld-wfs) < .01)
        break
    end
end

end

