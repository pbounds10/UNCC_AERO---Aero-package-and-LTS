function [ time, newGear ] = gearShift( currentGear, velo, gRatio, fd , wheelRad, rpmMax)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

time = 0;

fDrives = gRatio.*fd;

rpm = velo.*60.*fDrives./2./pi()./wheelRad;

for bb = 1:length(rpm)
    
    if (rpm < 3500)
        newGear = 1;
        break
    end 
    
    if (rpm(bb) < rpmMax(length(rpmMax)) && rpm(bb) > rpmMax(1) )
        newGear = bb;
        break
    end
end

if (currentGear ~= newGear)
    time = .2;
end




end

