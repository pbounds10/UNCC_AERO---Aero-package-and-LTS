function [ pos, dist ] = smoothTrack( pos, turn, maxDist, minDist, fitRange)

%This will add points to the track to ensure that atleast XX points can be
%made in every curve and that no point is more than YY meters way from each
%other
%   Detailed explanation goes here
% 
% pos = agumentTrackPos(pos); %This adds slight random numbers to the end of each x location in order to make sure that the curve fit is able to fit a curve to every location. If the exact same x location is present the fit will encounter errors
% 
% fitRange = fitPoints; %take 4 points ahead and behind of the given point
% 
% %calculate the initial distances
% for ii = 1:length(pos)
%     
%     if (ii < length(pos)) %if we haven't hit the last point then keep on going
%                 dist(ii) = pdist( [pos(ii,:); pos(ii+1,:)] ,'euclidean');
%     end 
%     
% end
% 
% %Start the adding of points
% posOld = [];
% while(1) %this will make the entire track iterate unitl every section satisifies the condition for minimum distance
%     
% %     if (isequal(posOld,pos))
% %         disp('error')
% %         pause
% %     end
%     posOld = pos;
%     distOld = dist;   
% for ii = 1:length(distOld) % now we will walk through the distance array and check each point to see if its exceeds the max distance
%     x = []; %clear the x and y data set
%     y = [];
%     bb = ii;
%     if (distOld(bb)-maxDist > .5) %This checks to see if the distance is to large
%        if (bb > fitRange && bb < length(posOld)-fitRange) %This checks to see if there are enough points ahead and behind
%           % disp('here1')
%            x = posOld((bb-fitRange):(bb+fitRange),1);
%            y = posOld((bb-fitRange):(bb+fitRange),2);
%            
%         elseif(bb <= fitRange && length(posOld)-fitRange)  %if the matrix is so short than the range is longer than the matrix then just use the whole set of points   
%           % disp('here4')
%            x = posOld(:,1);
%            y = posOld(:,2);
%           
%        elseif(bb <= fitRange) %Check to see if we are able to make points behind of the current
%            %disp('here2')
%            x = posOld((1):(bb+fitRange),1);
%            y = posOld((1):(bb+fitRange),2);
%            
%        elseif(bb >= length(distOld)-fitRange)%Check to see if we are able to make points ahead of the current
%            %disp('here3')
%            x = posOld((bb-fitRange):(length(posOld)),1);
%            y = posOld((bb-fitRange):(length(posOld)),2);
% 
%            
%        end
%  
%        
%        %Now fit a curve to the given data set
%        [curveFit,goodness] = fit(x,y, 'poly3');
%        leg1 = ['smooth ', num2str(goodness.sse)];
%        goodness.sse
%        
% %        if (goodness.sse > 1e-3)
% %            [curveFit,goodness] = fit(x,y, 'poly2');
% %                leg1 = ['poly3 ', num2str(goodness.sse)];
% %            if (goodness.sse > 1e-3)
% %                [curveFit,goodness] = fit(x,y, 'poly2');
% %                leg1 = ['poly2 ', num2str(goodness.sse)];
% %                if (goodness.sse > 1e-3)
% %                [curveFit,goodness] = fit(x,y, 'poly1');
% %                leg1 = ['poly1 ', num2str(goodness.sse)];
% %                end
% %            end
% %        end
% 
%        %Now use that curve to fill the space inbetween
%            numPoints = ceil(abs((posOld(bb,1)-posOld(bb+1,1)))/maxDist);
%            spacing = (posOld(bb+1,1)-posOld(bb,1))/(numPoints+1);
% 
%            
%            xNew = posOld(bb,1);
%            posNew = [];
%            
%            for kk = 1:numPoints              
%                xNew = xNew+spacing;                  %( (pos(ii,1)+pos(ii+1,1)) ) /2; %take the x location half way inbetween the bad point and the point after it
%                yNew = curveFit(xNew);  %get the y location
% %                if (abs(yNew-pos(bb,2)) > 20)
% %                    disp(['y error at', num2str(bb)])
% %                    curveFit
% %                    goodness.sse
% %                    plot(curveFit,posOld(1:end,1),posOld(1:end,2))
% %                    pos(bb-5:bb+5)
% %                    drawnow
% %                    
% %                end
%                posTemp = horzcat(xNew,yNew);            %make it a single matrix
%                posNew = [posNew; posTemp];
%                 
%            end
%            index = find(posOld(bb,1) == pos(:,1) & posOld(bb,2) == pos(:,2));
%            pos = cat(1,pos(1:index,:),posNew,pos((index+1):end,:)); %add the new point that we created into the position matrix
%           
%            if(length(index) > 1)
%                disp('index error')
%                pause
%            end
%            for kk = 1:length(pos) %Calculate the new distance between all points   
%                 if (kk < length(pos)) %if we haven't hit the last point then keep on going
%                             dist(kk) = pdist( [pos(kk,:); pos(kk+1,:)] ,'euclidean');
%                 end
%            end
%            clf
%            figure(1)
%            hold on
%            title(leg1)
%           % plot(curveFit,x,y)
%            scatter(posOld(1:end,1),posOld(1:end,2))
%            %scatter(posNew(1:end,1),posNew(1:end,2))
%            scatter(posOld(bb,1),posOld(bb,2),20,'fillded')
%            %scatter(posNew(:,1),posNew(:,2))
%            drawnow
%            %pause(1)
%            disp(sum(dist-maxDist > .1))
%     end
%     
% 
% % figure(2)
% % scatter(pos(:,1),pos(:,2))
% % figure(3)
% % scatter(pos(1:bb,1),pos(1:bb,2),80,'filled','d')
% % length(pos)
% % disp(sum(dist>1))
% 
% end
% 
% disp('out of for loop')
% 
% disp(sum(dist-maxDist > .5))
% 
%     if(dist-maxDist < .5) %if all points in the entire position matrix are below the threshold then end the entire program
%         disp('time to end the code')
%         break
%     end
% 
% end
% 
% 
% disp('end of code')
% % dist
% % pos
% % scatter(pos(:,1),pos(:,2))
% 
% end
% 

pos = agumentTrackPos(pos); %This adds slight random numbers to the end of each x location in order to make sure that the curve fit is able to fit a curve to every location. If the exact same x location is present the fit will encounter errors

[ numTurns, turnPosStart, turnPosEnd ] = defineCorners(pos,turn);

breakPoints = cat( 2, 1, turnPosStart, turnPosEnd );
breakPoints = sort( breakPoints);


curveFits = [];
xLocations(1) = pos(1,1);


ii = 0;
startInd = 1;
sets = zeros(length(pos) ,2);

for ii = 2:length(breakPoints)
     tempPos = pos( breakPoints(ii-1):breakPoints(ii), :);
     tempPos( length(tempPos):length(sets), : ) = 0;
     if (ii == 2)
         sets = tempPos;
         length(tempPos)
     else
        sets = cat( 3, sets, tempPos );
     end
end


ii = 0;

while(1)
    ii = ii + 1;
    
    realLength = max( max( sum( sets(:,:,ii) ~= 0, 1 ) ) );
    trimSets = sets(1:realLength,:,ii);
    
    curve{ii} = fit( trimSets(:,1),trimSets(:,2), 'smoothingSpline' );
    f = curve{ii};
    figure(ii)
    plot( f,sets(:,1,ii),sets(:,2,ii) )
    
    if ( ii == size(sets,3) )
        break
    end
end

maxDist = 2;
minDist = 1;

xNewPos = [];
yNewPos = [];

for ii = 1:size(sets,3)
    realLength = max( max( sum( sets(:,:,ii) ~= 0, 1 ) ) );
    trimSets = sets(1:realLength,:,ii);
    
    fx = curve{ii};
    %plot(fx)
    
    xStart = trimSets(1,1);
    xEnd = trimSets(length(trimSets),1);
    
    xNew = xStart + .1;
    
    %The distance formula using input function
    %fDist = @(x,y) sqrt( (x-y).^2 + (fx(x)-fx(y)).^2);
    %fDistEnd = @(y) sqrt( (xEnd-y).^2 + (fx(xEnd)-fx(y)).^2);
    %fDistEnd = distCalc( xEnd, fx(xEnd), xNew, fx(xNew) );
    while(1)
        
        xNew = findDistFunc(xStart,xNew,fx, maxDist, minDist);
        yNew = fx(xNew);
        xNewPos( size( xNewPos )+1 ) = xNew;
        yNewPos( size( yNewPos )+1 ) = yNew;
        
        fDistEnd = distCalc( xEnd, fx(xEnd), xNew, fx(xNew) );
        
        if( fDistEnd < maxDist )
            xNewPos( size( xNewPos )+1 ) = xEnd;
            yNewPos( size( yNewPos )+1 ) = fx(xEnd);
            break
        end
        
        xStart = xNew;
        
    end
%     figure(1)
%     scatter(trimSets(:,1),trimSets(:,2))
%     figure(2)
%     scatter(xNewPos,yNewPos)
    
    
end

sets

end
