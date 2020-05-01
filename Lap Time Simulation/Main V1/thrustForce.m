function [ thrust ] = thrustForce( velo, gRatio, fd, torque, wheelRad, rpmMAX  )
%This calculates the thrust force from the wheels based on the gear
%   Detailed explanation goes here

fDrives = gRatio.*fd;

wheelMulti = fDrives./wheelRad; %this is now in 1/m units

maxSpeed = rpmMAX(length(rpmMAX)).*2.*pi().*wheelRad./60./fDrives;

rpm = velo.*60.*fDrives./2./pi()./wheelRad;

gear = 0;
gearIndex = 0;

minRpm = 0;

if (rpm > rpmMAX(length(rpmMAX)) )
    disp("Error Transmission cannont achieve this speed")
    return
end

for bb = 1:length(rpm)
    
    if (rpm < 3500)

        gear = wheelMulti(1);
        gearIndex = 1;
        minRpm = 1;
        break
    end 
    
    if (rpm(bb) < rpmMAX(length(rpmMAX)) && rpm(bb) > rpmMAX(1) )
        gear = wheelMulti(bb);
        gearIndex = bb;
        break
    end
end

if (minRpm == 0)
    rpm = round(rpm(gearIndex)/500)*500;
else
    rpm = 3500;
end

rpmIndex = find(rpmMAX==rpm);


thrust = gear*torque(rpmIndex);


end

