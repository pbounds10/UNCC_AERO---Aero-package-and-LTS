function [u, ay, alphaF,alphaR,steerAngle,fyR,fyF] = findLatAccelMax(tireCo,camberF,camberR,curveRad,distF,distR,wF,wR,drag,driveWheels)
%This code finds the maximum lateral acceleration for a given radius
%Inputs:
%   tireCo - table of tire coefficients
%   camberF/R camber angle in degrees
%   curveRad - meters
%   dist F/R dist from axle front/rear to cg
%   wF/R - weight front/rear - N
%   driveWheels - rwd = 1,fwd = 2,awd =3
%Outputs:
%   u - forward velocity at max lat accel (m/s)
%   ay = lateral acceleration (m/s^2)
%   alphaF = slip angle front (deg)
%   alphaR = slip angle rear (deg)


%% ----- logic -----
%Using the goverening equations we assume that ax = 0 at the apex where
%this will be used. Then we assume find the peak tire force for the front
%and rear. then check to see which is the limiting force. the non-limiting
%force is then constrained in order to create a moment balence(yaw
%balence). then the force balence for the Y direction is solved to find u
%and ay. 
%To improve on this model you would assume that the rear tire produce
%enough force to overcome drag at this point. You would then back out the
%Kappa required to do this and use it to know how much grip has been taken
%away from the lateral grip of the tire because of the combined slip. (This
%does require a combined slip model for both Fx and Fy.

%% ----Initialize values
mass = (wF+wR)./9.81; %find mass in (kg)

%put weight in kN
wF = wF./1000/2; %per tire
wR = wR./1000/2; %per tire

%% -------find what force the drive wheels need to make to balence drag

if driveWheels == 1
    ft = @(kappa) 2.*magicTire(wR,0,kappa,camberR,tireCo,2) - drag;
    kNeg = .1;
    kPos = .1;
    while(ft(kNeg)>0)
        kNeg = kNeg - .1;
    end
    while(ft(kPos)<0)
        kPos = kPos + .1;
    end
    kappa = methodBisection(ft,kPos,kNeg);
end

%% -----Find peak lateral forces
%Get tire cornering stiffness and peak alpha values
    [~,~,aPeakF,~] = magicTire(wF,0,0,camberF,tireCo,1); %cornering stiffness is in N/deg
if driveWheels == 1
    [~,~,aPeakR,kappaPeak] = magicTire(wR,0,kappa,camberR,tireCo,1); %cornering stiffness is in N/deg
else
    [~,~,aPeakR,kappaPeak] = magicTire(wR,0,0,camberR,tireCo,1); %cornering stiffness is in N/deg
end

    
%Find peak lateral forces
    fyFPeak = 2.*magicTire(wF,aPeakF,0,camberF,tireCo,1);
if driveWheels == 1
    fyRPeak = 2.*magicTire(wR,aPeakR,kappa,camberR,tireCo,1);
else
    fyRPeak = 2.*magicTire(wR,aPeakR,0,camberR,tireCo,1);  
end

%% Find which lateral force is the limiting force
%litting force is the one which when mutliplied by its distance from cg is
%smaller

alpha_Array = 0:.01:max([aPeakF,aPeakR]);

if ( abs(fyFPeak.*distF) > abs(fyRPeak.*distR) )%this means that rear is limiting
    fyF = (distR./distF).*fyRPeak;
    fy_Array = 2.*magicTire(wF,alpha_Array,0,camberF,tireCo,1);
    
    alphaF = alpha_Array( find( min(abs(fy_Array - fyF))==(abs(fy_Array - fyF)) ) );
    alphaR = aPeakR;
    fyR = fyRPeak;
else   %this means that the front is the limiting
    fyR = (distF./distR).*fyFPeak;
    fy_Array = 2.*magicTire(wR,alpha_Array,0,camberR,tireCo,1);
    
    alphaR = alpha_Array( find( min(abs(fy_Array - fyR))==(abs(fy_Array - fyR)) ) );
    alphaF = aPeakF;
    fyF = fyFPeak;
end

%% solve for peak lateral acceleration
u = sqrt( ( curveRad.*(fyF+fyR) )./mass ); %forward velocity for max lat acceleration
ay = u.^2./curveRad; %lateral acceleration
steerAngle = alphaF + (distR./curveRad-alphaR) + distF./curveRad;

end

