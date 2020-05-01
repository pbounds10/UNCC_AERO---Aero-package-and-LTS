function [ thrust ] = thrustForce1( velo, gRatio, fd, torque, wheelRad, rpmArray  )
%This calculates the thrust force from the wheels based on the gear
%Inputs:
%velo - velocity m/s
%gRatio - array of gearRatios
%fd - final drive
%torque - array of torques that are from dyno data
%wheelRad - radius of drive wheels
%rpmArray - rpms at which dyno torque occured
%Outputs:
%thrust - force at wheels from engine N

fDrives = gRatio.*fd;

wheelMulti = fDrives./wheelRad; %this is now in 1/m units

maxSpeed = rpmArray(length(rpmArray)).*2.*pi().*wheelRad./60./fDrives;

rpm = velo.*60.*fDrives./2./pi()./wheelRad;

gear = 0;
gearIndex = 0;

minRpm = 0;

% if (rpm > rpmMAX(length(rpmMAX)) )
%     rpm = rpmMAX(length(rpmMAX));
%     disp("Error Transmission cannont achieve this speed. Clipping to max RPM")
% end

for bb = 1:length(rpm)
    
    if (rpm < rpmArray(1)) %if all rpms are below the min then set gear = 1;

        gear = wheelMulti(1);
        gearIndex = 1;
        minRpm = 1;
        break
    end 
    
    if (rpm(bb) < rpmArray(length(rpmArray)) && rpm(bb) > rpmArray(1) ) %check for the 1st gear that is less than maxRpm but greater than minRpm 
        gear = wheelMulti(bb);
        gearIndex = bb;
        break
    end
end

if (minRpm == 0) %if the condition of minRpm is not true then just calculate like normal
    rpm = round(rpm(gearIndex)/500)*500; %round to the nearest 500 rpm
else %else set to min
    rpm = rpmArray(1);
end

rpmIndex = find(rpmArray==rpm);


thrust = gear*torque(rpmIndex);


end

