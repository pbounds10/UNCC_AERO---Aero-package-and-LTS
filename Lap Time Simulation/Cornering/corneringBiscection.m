function [uMax,ay,alphaF,alphaR,steerAngle,driveGrip] = corneringBiscection(tireCo,camberF,camberR,curveRad,distF,distR,wF,wR,maxSteer,driveWheels)
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

%% ----Initialize values ---- %%
mass = (wF+wR)./9.81; %find mass in (kg)
optTol = .01;

%put weight in kN
wF = wF./1000/2; %per tire
wR = wR./1000/2; %per tire


alphaF_store = [];
alphaR_store = [];

kk = 1; %array counter

%uu = 1:.5:70; %make a velocity matrix for sweeping
%uu = uu';

steerAngleArray = 0:.25:maxSteer; %make a steerAngle array to chose the max value from

%% ---- Start steering angle sweep ---%%%
for steerAngle = steerAngleArray
    alphaF = [];
    alphaR = [];

        [alphaF,alphaR] = yawSettle(steerAngle,wF,wR,curveRad,tireCo,distR,distF,camberF,camberR); 
         %If a negative alphaF value occurs this means that the steering angle
         %was not enough to complete the turn aka it was lower than the akerman
         %angle.
            if (isnan(alphaF))
                continue
            end
       
        fYF = 2.*magicTire(wF,alphaF,0,camberF,tireCo,1); %cornering stiffness is in N/deg
        fYR = 2.*magicTire(wR,alphaR,0,camberR,tireCo,1); %cornering stiffness is in N/deg
        
        vD = @(uu) ( abs(( fYF + fYR )./mass) - (uu.^2)./curveRad ); %find lateral accerlation
        
        if(vD(0) < 0 ) %check to see if there are any solutions
            continue
        end
        
        alphaF_store(kk) = alphaF;
        alphaR_store(kk) = alphaR;
        fyR_store(kk) = fYR;
        fYF_store(kk) = fYF;
        steadyCheck(kk) = fYF.*distF-fYR.*distR;
        
        
        %% ---- start bisection method
            uNeg = 1;
            uPos = .1;
            yNeg = vD(uNeg); %make a value to be negative
            yPos = vD(uPos); %make a very small value that is sure to be positive    
            
            while(sign(yNeg) == 1)
                 yNeg = vD(uNeg);
                 uNeg = uNeg + 5;
            end
            
            uNew = 1;
        
         while(1)
             xOld = uNew;

             uNew = (uPos+uNeg)./2;
             
             if (vD(uNew) < 0 )
                 uNeg = uNew;
             else
                 uPos = uNew;
             end
            
             if (abs(uNew-xOld) < .00001)
                 break
             end
             
         end     
             
        uMax(kk) = uNew;
        
        ay(kk) = uMax(kk).^2./curveRad;
        
        kk = kk + 1;

end

%% --- find the peak values
indMax = find(ay == max(ay),1,'first'); %find the position of the max lat Accel

ay = ay(indMax);
uMax = uMax(indMax);
steerAngle = steerAngleArray(indMax);
alphaF = alphaF_store(indMax);
alphaR = alphaR_store(indMax);



%% ---- Calculate grip left to accelerate in the long direction
if (driveWheels == 1)
    %Get peak kappa values
    [~,~,~,kPeakR] = magicTire(wR,alphaR,0,camberR,tireCo,2);

    %Find peak lateral forces
    driveGrip = 2.*magicTire(wR,alphaR,kPeakR,camberR,tireCo,2); %fx force avaliable in newtons
    
end




end

