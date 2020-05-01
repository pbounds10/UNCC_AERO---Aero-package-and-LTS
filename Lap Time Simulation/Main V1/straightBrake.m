function [ posInd, err1 ] = straightBrake( vTrack, vCurv, cStart, cEnd) %dist, accel, curentPos )
%this checks to see if there is a place that merges straights and curves
%   takes a velocity profile from the straight and curve
%   takes the index location of the start and end of the curve


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DOES NOT NEED TO BE CHANGED NO PHYSICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALCULATED IN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lStr = length(vStr);
% lCurv = length(vCurv);
err1 = 0;

tol = .25;

    
% for ii = cStart:-1:cEnd %if we havent hit the end of the track do this
%            if (ii == cStart) %if we are at the start position use the max Turn velocity
%                 vCurv(ii) = constABraking(accel, vTrack(ii), dist(ii));
%            else
%                 vCurv(ii) = constABraking(accel, vCurv(ii+1), dist(ii)); %it's plus 1 because we are going backwards  
%            end
% end 

for ii = (cStart - 1):-1:cEnd
    
    
    
    if( abs(vTrack(ii) - vCurv(ii)) <= tol)
        posInd = ii;
         %disp(num2str(abs(vTrack(ii) - vCurv(ii))) + " right on braking speed" );
        return
    elseif( vCurv(ii) > vTrack(ii) )
        posInd = ii;
        %disp(num2str(vCurv(ii)) + " braking");
        %disp(vTrack(ii));
        return
    end
    
end

err1 = 1; % this tells the system that you could not find a point that let you merge the two
posInd = NaN;




end

