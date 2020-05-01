function [uMax,ay,alphaF,alphaR,steerAngle,driveGrip] = cornering(tireCo,camberF,camberR,curveRad,distF,distR,wF,wR,maxSteer,driveWheels)
%This code finds the lateral and long accel in a corner for a given radius
%Inputs:
%   tireCo - table of tire coefficients
%   camberF/R camber angle in degrees
%   curveRad - meters
%   dist F/R dist from axle front/rear to cg
%   wF/R - weight front/rear - N
%   SteerMax - max lock on wheels (deg)
%   driveWheels - rwd = 1,fwd = 2,awd =3
%Outputs:
%   u - forward velocity at max lat accel (m/s)
%   ay = lateral acceleration (m/s^2)
%   alphaF = slip angle front (deg)
%   alphaR = slip angle rear (deg)

%% ----Initialize values
mass = (wF+wR)./9.81; %find mass in (kg)
optTol = .01;

%put weight in kN
wF = wF./1000/2; %per tire
wR = wR./1000/2; %per tire


alphaF_store = [];
alphaR_store = [];

kk = 1; %array counter

uu = 1:.5:70; %make a velocity matrix for sweeping
uu = uu';

steerAngleArray = 0:.25:maxSteer; %make a steerAngle array to chose the max value from

%for the given curve do a steering angle sweep in deg
for steerAngle = steerAngleArray
    alphaF = [];
    alphaR = [];

        [alphaF,alphaR] = yawSettle(steerAngle,wF,wR,curveRad,tireCo,distR,distF,camberF,camberR);              
       
        fYF = 2.*magicTire(wF,alphaF,0,camberF,tireCo,1); %cornering stiffness is in N/deg
        fYR = 2.*magicTire(wR,alphaR,0,camberR,tireCo,1); %cornering stiffness is in N/deg
        
        vD(:,1) = ( (( fYF + fYR )./mass) - (uu.^2)./curveRad ); %find lateral accerlation
        
        if(sum( sign(vD) == -1) == size(vD,1) )
            %disp("Steering angle too low, no possible solutions")
            continue
        end
       
        alphaF_store(kk) = alphaF;
        alphaR_store(kk) = alphaR;
       
        fVD = fit(uu,vD,'poly3');
        
        %find the numerical differential of vD
        for ii = 1:length(uu)
            if (ii < length(uu)-2)
                vDD(ii,1) = unequalDiff(uu(ii:ii+2),fVD(uu(ii:ii+2)),uu(ii));
            else
                vDD(ii,1) = unequalDiff(uu(ii-2:ii),fVD(uu(ii-2:ii)),uu(ii));
            end
        end
                
        fVDD = fit(uu,vDD,'poly1'); %fit a curve to the data
        
        uMax(kk) = newtonRaphson(fVD,fVDD,uu(1),optTol);
        
        ay(kk) = uMax(kk).^2./curveRad;
        kk = kk + 1;

end


indMax = find(ay == max(ay),1,'first'); %find the position of the max lat Accel

ay = ay(indMax);
uMax = uMax(indMax);
steerAngle = steerAngleArray(indMax);
alphaF = alphaF_store(indMax);
alphaR = alphaR_store(indMax);

if (driveWheels == 1)
    %Get peak kappa values
    [~,~,~,kPeakR] = magicTire(wR,alphaR,0,camberR,tireCo,2);

    %Find peak lateral forces
    driveGrip = 2.*magicTire(wR,alphaR,kPeakR,camberR,tireCo,2); %fx force avaliable in newtons
    
end




end

