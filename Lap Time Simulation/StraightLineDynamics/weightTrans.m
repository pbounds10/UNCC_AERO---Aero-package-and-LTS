function [ fFront, fRear ] = weightTrans( wfs, wrs, cgH, len,lenA, accel )
%calculates the longitudinal weight transfer due to acceleration
%Inputs:
%   wfs - static weight of front
%   wrs - static weight of rear
%   cgH - height of the center of gravity in m
%   len - length of the car in m
%   accel - acceleration of car m/s^2
%   lenA - distance from front axle to cg
%Output:
%   fFront - force on front (N)
%   fRear - force on rear (N)
g = 9.81; %   g - acceleration due to gravity m/s^2

mass = (wfs./g+wrs./g); % total mass in kg

lenB = len - lenA; %length to rear axle in m

fFront = (lenB.*mass.*g)./(len)-(cgH.*mass.*accel)./(len);
fRear = (lenA.*mass.*g)./(len)+(cgH.*mass.*accel)./(len);

%since springs are not included to stop full weight transfer put in a max
%transfer
% if (fFront > 0.8*(wfs+wrs))
%     fFront = 0.8.*(wfs+wrs);
%     fRear = 0.2.*(wfs+wrs);
% elseif(fRear > 0.8*(wfs+wrs))
%     fFront = 0.2.*(wfs+wrs);
%     fRear = 0.8.*(wfs+wrs);
% end


end

