function [ax,wfs,wrs,kappaPeak] = brakingAccel(u,v, pressIn, pressMax, radWheelF,radWheelR, bgf, bgr, pCrit, pro, wfs, wrs, cgH, len, lenA, SA, tireCos, camberF, camberR, springF,springR, yaw, pitch, aeroCD, aeroCDF, aeroFront, area, density, symetric, curveRad)
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


%% Loop til steady
while(1) 
    %% 1) find max traction for braking
        [~,~,~,kappaPeak] = magicTire(wrs./2000,0,0,camberR,tireCos,2);
        tireForceRear = 2.*magicTire(wrs./2000,SA,kappaPeak,camberR,tireCos,2);
        [~,~,~,kappaPeak] = magicTire(wfs/2000,0,0,camberF,tireCos,2);
        tireForceFront = 2.*magicTire(wfs/2000,SA,kappaPeak,camberF,tireCos,2);      
    
    %% 2) find max braking force from brake system
        [fFront,fRear] = brakeByPress( pressIn, pressMax, wfs, wrs, radWheelF, radWheelR, bgf, bgr, pCrit, pro);
        
    %% 3) find limit brakes or tires
    if (fFront > tireForceFront)
        FxFront = tireForceFront;
    else
        FxFront = fFront;
    end
    if (fRear > tireForceRear)
        FxRear = tireForceRear;
    else
        FxRear = fRear;
    end
    
    Fx = FxFront +FxRear;

    %% 4) find acceleration
    drag = aeroForce(yaw,pitch,u,v,aeroCD,aeroCDF,aeroFront,area,density,1);
    if (~isnan(curveRad))
        ax = (-Fx-drag)./mass -(v.*u)./curveRad;
    else
        ax = (-Fx-drag)./mass;
    end

    %% 5) find weight transfer
    wfsOld = wfs;
    [wfs,wrs] = weightTrans(wfs, wrs, cgH, len,lenA, ax );
    
    %% 6) find ride height change and pitch for aero
        %currently not do but could be

    %% 7) check to see if its time to break loop
    if (abs(wfsOld-wfs) < .01)
        break
    end
end

end

