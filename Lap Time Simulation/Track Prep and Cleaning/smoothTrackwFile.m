function [ pos, dist ] = smoothTrackwFile( segmentFile, maxDist, minDist, fitRange)
%This will add points to the track to ensure that atleast XX points can be

[x, y, segment] = importfileSegment(segmentFile, "Sheet1", [2, 1000]);

x(isnan(x) ) = [];
y(isnan(y) ) = [];
segment(isnan(segment) ) = [];

pos = horzcat(x,y);

pos = agumentTrackPos(pos); %This adds slight random numbers to the end of each x location in order to make sure that the curve fit is able to fit a curve to every location. If the exact same x location is present the fit will encounter errors

pointRemove = [];

for ii = 2:size(pos,1)-1
    if ( pos(ii-1,1) > pos(ii,1) && pos(ii+1,1) > pos(ii,1) )
        pointRemove = cat(1,pointRemove,ii);
        disp('removed point due to min error')
    elseif ( pos(ii-1,1) < pos(ii,1) && pos(ii+1,1) < pos(ii,1) )
        pointRemove = cat(1,pointRemove,ii);
        disp('removed point due to max error')
    end
    
    if ( pos(ii-1,2) > pos(ii,2) && pos(ii+1,2) > pos(ii,2) )
        pointRemove = cat(1,pointRemove,ii);
        disp('removed point due to min error')
    elseif ( pos(ii-1,2) < pos(ii,2) && pos(ii+1,2) < pos(ii,2) )
        pointRemove = cat(1,pointRemove,ii);
        disp('removed point due to max error')
    end
end

pos(pointRemove,:) = [];
segment(pointRemove) = [];




[ numTurns, turnPosStart, turnPosEnd ] = defineCorners(pos,segment);

breakPoints = cat( 2, 1, turnPosStart, turnPosEnd );
breakPoints = sort( breakPoints);


curveFits = [];
xLocations(1) = pos(1,1);


ii = 0;
startInd = 1;
sets = {};%zeros(length(pos) ,2);

for ii = 2:length(breakPoints)
     tempPos = {pos( breakPoints(ii-1):breakPoints(ii), :)};
     %tempPos( length(tempPos):length(sets), : ) = 0;

        sets = cat( 1, sets, tempPos );

end


ii = 0;

while(1)
    ii = ii + 1;
    
    %realLength = max( max( sum( sets(:,:,ii) ~= 0, 1 ) ) );
    trimSets = cell2mat( sets(ii) ); %sets(1:realLength,:,ii);
    curve{ii} = fit( trimSets(:,1),trimSets(:,2), 'smoothingSpline' );
%     f = curve{ii};
%     figure(ii)
%     plot( f,sets(:,1,ii),sets(:,2,ii) )
    
    if ( ii == size(sets,1) )
        break
    end
end


xNewPos = [];
yNewPos = [];

for ii = 1:size(sets,1)
    %realLength = max( max( sum( sets(:,:,ii) ~= 0, 1 ) ) );
    trimSets = cell2mat( sets(ii) );
    %inc = mean( diff(abs(trimSets(:,1))) );
    trimSets(3:size(trimSets,1)-1,:) = [];
    fx = curve{ii};
    %plot(fx)
    
    
    xStart = trimSets(1,1);
    yStart = trimSets(1,2);
    disp(ii)
    xEnd = trimSets(length(trimSets),1);
    
    signAugX = sign(xEnd-xStart);

        inc = roundsd((xEnd - xStart),1)/10;

    xNext = xStart+signAugX.*abs(inc); %xStart + (trimSets(2,1) - xStart)./2;
    %inc = .1;
    
    %The distance formula using input function
    %fDist = @(x,y) sqrt( (x-y).^2 + (fx(x)-fx(y)).^2);
    %fDistEnd = @(y) sqrt( (xEnd-y).^2 + (fx(xEnd)-fx(y)).^2);
    %fDistEnd = distCalc( xEnd, fx(xEnd), xNew, fx(xNew) );
    
    zz = 1;
    nextPos = [];
    xNew = [];
    yNew= [];
    tempX = [xStart];
    
    while(1)
        [xNew,kk] = findDistFunc(xStart,xNext,fx, maxDist, minDist, inc);
        yNew = fx(xNew);
        xNewPos( size( xNewPos )+1 ) = xNew;
        yNewPos( size( yNewPos )+1 ) = yNew;
        
        tempX( size( tempX )+1 ) = xNew;
        
        fDistEnd = distCalc( xEnd, fx(xEnd), xNew, fx(xNew) )
        
        if( fDistEnd < maxDist )
            xNewPos( size( xNewPos )+1 ) = xEnd;
            yNewPos( size( yNewPos )+1 ) = fx(xEnd);
            break
        end
        
%         signAugX = sign(xNew-xStart);
%         signAugY = sign(yNew-yStart);
%         if signAugX == 1
%             if( sum( sign( ( trimSets(:,1)-xNew ) ) ) == size(trimSets,1) )
%                 if ( signAugY == 1 )
%                     nextPos(zz) = find( sign( ( trimSets(:,2)-yNew ) ) == 1,1)
%                 else
%                 	nextPos(zz) = find( sign( ( trimSets(:,2)-yNew ) ) == -1,1)
%                 end
%             end
%             nextPos(zz) = find( sign( ( trimSets(:,1)-xNew ) ) == 1,1)
%         elseif signAugX == -1
%             
%             if(sum( sign( ( trimSets(:,1)-xNew ) ) ) == size(trimSets,1) )
%                 if ( signAugY == 1 )
%                     nextPos(zz) = find( sign( ( trimSets(:,2)-yNew ) ) == 1,1)
%                 else
%                 	nextPos(zz) = find( sign( ( trimSets(:,2)-yNew ) ) == -1,1)
%                 end
%             end
% 
%             nextPos(zz) = find( sign( ( trimSets(:,1)-xNew ) ) == -1,1)
%          end
        inc = roundsd((xEnd - xNew),1)

        xStart = xNew;
        yStart = yNew;
        xNext = xStart+signAugX.*abs(inc); %trimSets(zz,1)
        zz = zz + 1;
    end
    figure(1)
    scatter(trimSets(:,1),trimSets(:,2))
    figure(2)
    scatter(xNewPos,yNewPos)
    
    
end

sets

end
