function [alphaF,alphaR] = yawSettle(steerAngle,wF,wR,curveRad,tireCo,distR,distF,camberF,camberR)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

kk = 1;

alphaR = 0; %find the rear slip angle in deg
sideSlip = 57.3.*distR./curveRad - alphaR; %find the side slipe in deg
alphaF = steerAngle - sideSlip - 57.3.*distF./curveRad ;
error = 1000;
while(1)
    error_old = error;
    alphaR_old = alphaR;
    alphaF_old = alphaF;
    
    sideSlip = 57.3.*distR./curveRad - alphaR; %find the side slipe in deg
    alphaF = steerAngle - sideSlip - 57.3.*distF./curveRad ;
    %If a negative alphaF value occurs this means that the steering angle
    %was not enough to complete the turn aka it was lower than the akerman
    %angle.
    if(alphaF < 0)
        alphaF = NaN;
        alphaR = NaN;
        return
    end

    fYF = 2.*magicTire(wF,alphaF,0,camberF,tireCo,1); %cornering stiffness is in N/deg

    fYR = (distF./distR).*fYF; %putting in the constraint that we are in steady state conering yaw accel = 0


    alpha_sweep = -5:.1:20;
    
    %Find the peak fYR and see if fYF may have to be pull back
    [~,~,alphaPeakR] = magicTire(wR,0,0,camberR,tireCo,1);
    fYR_Peak = 2.*magicTire(wR,alphaPeakR,0,camberR,tireCo,1);
    
    if (fYR > fYR_Peak)
       fYF = (distR./distF).*fYR_Peak;
       fYR = fYR_Peak;
       
       %find the alpha value for the front
       fYF_curve = 2.*magicTire(wF,alpha_sweep,0,camberF,tireCo,1);
       fYF_diff = fYF_curve - fYF;
       indexF = find(abs(fYF_diff) == min(abs(fYF_diff)));
       
       alphaF = alpha_sweep( indexF );
       alphaR = alphaPeakR;
       
       error = ( distF.*fYF_curve(indexF) ) - ( distR.*fYR );
    else
        %Find the alpha value for the rear
        fYR_curve = 2.*magicTire(wR,alpha_sweep,0,camberR,tireCo,1);
        fYR_diff = fYR_curve - fYR;
        indexR = find(abs(fYR_diff) == min(abs(fYR_diff)));

        alphaR = alpha_sweep( indexR );
        error = ( distF.*fYF ) - ( distR.*fYR_curve( indexR ) );
    end

    
    
    if(abs(error) < abs(error_old))
        break
    end
    
    kk = kk + 1;
    
end

end

