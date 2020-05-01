function [ diffs, turn, slopes ] = characterizeTrack( pos, curveLen, slopeTol )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

 diffs = zeros( length(pos), 2);
% 
% for ii = 1:length(pos)
% 
%     if( ii == 1)
%     xDiff = pos(ii:ii+2,1);
%     yDiff = pos(ii:ii+2,2);
%     diffs(ii,1) = unequalDiff(xDiff,yDiff,xDiff(1));
%     diffs(ii,2) = unequalDiff(yDiff,xDiff,yDiff(1));
%     elseif(ii == length(pos))
%     xDiff = pos(ii-2:ii,1);
%     yDiff = pos(ii-2:ii,2);
%     diffs(ii,1) = unequalDiff(xDiff,yDiff,xDiff(3));
%     diffs(ii,2) = unequalDiff(yDiff,xDiff,yDiff(3));
%     else
%     xDiff = pos(ii-1:ii+1,1);
%     yDiff = pos(ii-1:ii+1,2);
%     diffs(ii,1) = unequalDiff(xDiff,yDiff,xDiff(2));
%     diffs(ii,2) = unequalDiff(yDiff,xDiff,yDiff(2));
%     end
%     
% end
% 
% diffs = abs( diffs );

% xIndAll = find( diffs(:,1) == Inf | diffs(:,1) == -Inf );
% yIndAll = find( diffs(:,2) == Inf | diffs(:,2) == -Inf );
% 
% xIndAll = cat( 1,xIndAll,1 );
% yIndAll = cat( 1,yIndAll,1 );
% 
% xIndLogic = diff(xIndAll) == 1;
% yIndLogic = diff(yIndAll) == 1;
% 
% xNum = 1;
% 
% for ii = 1:length( xIndLogic ) 
%     
%     if (ii < length(xIndLogic) )
%         if ( xIndLogic(ii) == 1 && xIndLogic(ii+1) == 0)
%             xInd(xNum) = xIndAll(ii);
%             xNum = xNum + 1;
%         end
%     end
%     
% end
     
turn = zeros( length( pos ),1 );

% for ii = xInd
%    turn( ii-curveLen:ii) = 1;   
%     
% end


counter = 1;

fitRange = 3;

for ii = 1:length(pos)
    
    if (ii <= fitRange)
        curve1 = fit( pos(ii:ii+fitRange,1), pos(ii:ii+fitRange,2), 'poly1');
    elseif (ii >= length(pos)-fitRange)
        curve1 = fit( pos(ii-fitRange:length(pos),1), pos(ii-fitRange:length(pos),2), 'poly1');
    else
        curve1 = fit( pos(ii-fitRange:ii+fitRange,1), pos(ii-fitRange:ii+fitRange,2), 'poly1');
    end
    
%     figure(1)
%     plot(curve1, pos(ii-3:ii+3,1), pos(ii-3:ii+3,2))
%     figure(2)
%     scatter(pos(ii-3:ii+3,1),pos(ii-3:ii+3,2) )
%     drawnow
%     pause(.5)

    slopes(counter,1) = curve1.p1;
    counter = counter + 1;
    
end

slopes( abs( slopes) > 1000 ) = 1000.*sign( slopes ( logical( abs(slopes) > 1000 )));

inflectionInd = findInflection(slopes);

xInd = find( diff ( sign (slopes) ) );

xInd = unique( sort( cat( 1, xInd, inflectionInd ) ) );

xInd( abs( diff(xInd) ) == 1 ) = [];

xInd( abs( slopes(xInd) ) < 1 ) = [];

turnGroup = findCurveStart(xInd, slopes)';

turn = zeros( length( pos ),1 );

for ii = 1:length(turnGroup)
    
    turn( turnGroup(ii) ) = 1;
    
end





% straightGroup = find ( abs( diff ( slopes ) ) < slopeTol);
% 
% turnGroup = find ( abs( diff ( slopes ) ) > slopeTol);
% 
% straightLog = diff( find ( abs( diff ( slopes ) ) < slopeTol) ) == 1;
%     
% while(1)
%     temp = find (straightLog == 0,2);
%     
%     if( diff( temp ) < 5)
%         straightGroup( temp(1):temp(2) ) = [];  
%     end
%         
%     if ( isempty( temp ) )
%         break
%     else
%         straightLog( temp(1):temp(2) ) = [];
%     end
%     
% end
% 
% 
% 
% xInd( ismember(xInd, straightGroup) == 1 ) = [];
% 
% turn = zeros( length( pos ),1 );
% 
% for ii = 1:length(xInd)
%     
%     indexStart = find( turnGroup == xInd(ii) );
%     
%     temp = logical( diff ( turnGroup(1:indexStart) ) < 3 );
%     indexEnd = find( temp == 0, 1, 'last' );
%     
%     indexStart = turnGroup(indexStart);
%     indexEnd = turnGroup(indexEnd);
%     
%     turn(indexEnd:indexStart) = 1;
% 
% 
% end    
    


    



end

