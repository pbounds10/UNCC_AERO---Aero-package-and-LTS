function [ turnPosStart, apexStart, apexEnd, rads, apexPos ] = defineCorners( pos, turn )
%This defines each turn on the track to get its parameters
%
%Logic: 
%Sort the track to find what places have been considered turns. place them
%in their own array and then start find their start point, end points, and
%rad
turnCol = turn;
turn = find(turn == 1);
divs = [0 find(diff(turn)~=1)']; %find the divisions in the data

numTurns = length(divs);
apexPos = [];
apexStart = [];
apexEnd = [];

rads = NaN(size(pos,1),1); %if it is a straight then no radius

kk = 1;

t = {};

for ii = 1:length(divs)
    if (ii ~= length(divs))
        t{kk} = {turn(divs(ii)+1:divs(ii+1))};
        turnPosStart(kk) = turn(divs(ii)+1);
        if (size(t{1,kk}{1},1) < 8)
            continue
        end

        turnPosEnd(kk) = turn(divs(ii+1));
        set = pos(t{1,kk}{1,1},:);
        circ = CircleFitByPratt(set);
        %rads(kk) = circ(3); %get the radius of the turn
    else
        t{kk} = {turn(divs(ii)+1:length(turn))};
        turnPosStart(kk) = turn(divs(ii)+1);
        if (size(t{1,kk}{1},1) < 6)
            continue
        end

        turnPosEnd(kk) = turn(length(turn));
        
        set = pos(t{1,kk}{1,1},:);
        
        circ = CircleFitByPratt(set);
        %rads(kk) = circ(3); %get the radius of the turn
    end
    kk = kk + 1;
end

%% section up each turn into several segments and assign radius to the segments

segmentLen = 10; %The optimal number of points in a segment
divArray = segmentLen-3:segmentLen+5; %Check if no optimal spacing is needed
indStart = 1;
numSets = 1;

for ii = 1:size(t,2) %for each set
    setAll = pos(t{1,ii}{1,1},:);
    indexs = t{1,ii}{1,1};
    segmentLen = 10;
    numSets = 1;
    indStart = 1;
    while(round(numSets) <=3)
        numSets = size(t{1,ii}{1,1},1)./segmentLen;
        setLen = segmentLen;
        divArry = segmentLen-3:segmentLen+5;
        if ((numSets-floor(numSets)).*segmentLen < 3) %check to see if a segment has less than 5 points
            numSets = size(t{1,ii}{1,1},1)./divArray;
            setLen = divArray(find((numSets-floor(numSets)).*segmentLen>5,1));
            numSets = numSets(find((numSets-floor(numSets)).*segmentLen>5,1));       
        end
        segmentLen = segmentLen - 1;
    end
    %% assign radii
	for aa = 1:round(numSets)
        if (aa~=round(numSets))
            tempSet = setAll(indStart:setLen.*aa,:);
            circ = CircleFitByPratt(tempSet);
            tempRad = circ(3);
            setLenMod = setLen.*aa + 3;
            setStartMod = indStart-3;
            while(tempRad > 400 || isnan(tempRad))   
               if (setLenMod < size(setAll,1))
                   tempSet = setAll(indStart:setLenMod,:);
                   circ = CircleFitByPratt(tempSet);
                   tempRad = circ(3);
                   setLenMod = setLenMod + 3;
               else
                   setLenMod = size(setAll,1);
                   tempSet = setAll(setStartMod:setLenMod,:);
                   circ = CircleFitByPratt(tempSet);
                   tempRad = circ(3);
                   setStartMod = setStartMod-3;
               end
            end
               
            rads(indexs(indStart:setLen.*aa)) = tempRad;
%             figure(1)
%             scatter(pos(:,1),pos(:,2),10,turnCol)
%             hold on
%             scatter(tempSet(:,1),tempSet(:,2))
        else
            tempSet = setAll(indStart:size(setAll,1),:);
            circ = CircleFitByPratt(tempSet);
            tempRad = circ(3);
            setStartMod = indStart-3;
            while(tempRad > 200 || isnan(tempRad))
               tempSet = setAll(setStartMod:size(setAll,1),:);
               circ = CircleFitByPratt(tempSet);
               tempRad = circ(3);
               setStartMod = setStartMod-3;
            end
            rads(indexs(indStart:size(indexs,1))) = tempRad;
%             tempRad
%             figure(1)
%             scatter(pos(:,1),pos(:,2),10,turnCol)
%             hold on
%             scatter(tempSet(:,1),tempSet(:,2))
        end
        indStart = setLen.*aa + 1;
        clf
    end
    %% check for high spots
    indStart = 1;
    for aa = 1:round(numSets)-2
         if (aa~=round(numSets)-2)
            if( mean(rads(indexs(indStart:setLen.*aa)) ) < mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))) & mean(rads(indexs(setLen.*(aa+1)+1:setLen.*(aa+2))) ) < mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))))%does the radius spikes and then return down low?
                tempRad = ( mean(rads(indexs(indStart:setLen.*aa)) ) + mean(rads(indexs(setLen.*(aa+1)+1:setLen.*(aa+2))) ) )./2;
                rads(indexs(setLen.*aa+1:setLen.*(aa+1))) = tempRad;
            end
              
         else
           if( mean(rads(indexs(indStart:setLen.*aa)) ) < mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))) & mean(rads(indexs(setLen.*(aa+1)+1:size(indexs,1))) ) < mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))))%does the radius spikes and then return down low?
                tempRad = ( mean(rads(indexs(indStart:setLen.*aa)) ) + mean(rads(indexs(setLen.*(aa+1)+1:size(indexs,1))) ) )./2;
                rads(indexs(setLen.*aa+1:setLen.*(aa+1))) = tempRad;

            end
        end
        indStart = setLen.*aa + 1;    
    end 
    %% check for low spots
    indStart = 1;
       for aa = 1:round(numSets)-2
         if (aa~=round(numSets)-2)
            if( mean(rads(indexs(indStart:setLen.*aa)) ) > mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))) & mean(rads(indexs(setLen.*(aa+1)+1:setLen.*(aa+2))) ) > mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))))%does the radius spikes and then return down low?
                tempRad = ( mean(rads(indexs(indStart:setLen.*aa)) ) + mean(rads(indexs(setLen.*(aa+1)+1:setLen.*(aa+2))) ) )./2;
                rads(indexs(setLen.*aa+1:setLen.*(aa+1))) = tempRad;
            end
              
         else
           if( mean(rads(indexs(indStart:setLen.*aa)) ) > mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))) & mean(rads(indexs(setLen.*(aa+1)+1:size(indexs,1))) ) > mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))))%does the radius spikes and then return down low?
                tempRad = ( mean(rads(indexs(indStart:setLen.*aa)) ) + mean(rads(indexs(setLen.*(aa+1)+1:size(indexs,1))) ) )./2;
                rads(indexs(setLen.*aa+1:setLen.*(aa+1))) = tempRad;

            end
        end
        indStart = setLen.*aa + 1;    
       end 
       
       %% check for large jumps
       indStart = 1;
       minRadJump = 10; %largest difference in sequential radii
       for aa = 1:round(numSets)-1
           if (aa == 1)
               if ( mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))) - mean(rads(indexs(indStart:setLen.*aa)) )< -minRadJump )
                   tempRad = mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))) + minRadJump;
                   rads(indexs(indStart:setLen.*aa)) = tempRad;
               end
               
           elseif (aa~=round(numSets)-1 & aa~=1)
               if ( mean(rads(indexs(indStart:setLen.*aa)) ) -  mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))) < -minRadJump )
                   tempRad = mean(rads(indexs(indStart:setLen.*aa))) + minRadJump;
                   rads(indexs(setLen.*aa+1:setLen.*(aa+1))) = tempRad;
               end
           else
               if ( mean(rads(indexs(indStart:setLen.*aa)) ) -  mean(rads(indexs(setLen.*aa+1:size(indexs,1)))) < -minRadJump )
                   tempRad = mean(rads(indexs(indStart:setLen.*aa))) + minRadJump;
                   rads(indexs(setLen.*aa+1:size(indexs,1))) = tempRad;
               end
           end
           indStart = setLen.*aa + 1;   
       end
       
       %% check for large drops
       indStart = 1;
       minRadDrop = 15; %largest difference in sequential radii
       for aa = 1:round(numSets)-1
           if (aa == 1)
               if ( mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))) - mean(rads(indexs(indStart:setLen.*aa)) )> minRadDrop )
                   tempRad = mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))) - minRadJump;
                   rads(indexs(indStart:setLen.*aa)) = tempRad;
               end
               
           elseif (aa~=round(numSets)-1 & aa~=1)
               if ( mean(rads(indexs(indStart:setLen.*aa)) ) -  mean(rads(indexs(setLen.*aa+1:setLen.*(aa+1)))) > minRadDrop )
                   tempRad = mean(rads(indexs(indStart:setLen.*aa))) - minRadJump;
                   rads(indexs(setLen.*aa+1:setLen.*(aa+1))) = tempRad;
               end
           else
               if ( mean(rads(indexs(indStart:setLen.*aa)) ) -  mean(rads(indexs(setLen.*aa+1:size(indexs,1)))) > minRadDrop )
                   tempRad = mean(rads(indexs(indStart:setLen.*aa))) - minRadJump;
                   rads(indexs(setLen.*aa+1:size(indexs,1))) = tempRad;
               end
           end
           indStart = setLen.*aa + 1;   
       end
    %% find the apex set
        minRad = indexs(find( rads(indexs) == min(rads(indexs)) ));
        if (max(diff(minRad))>1)
            minRad = minRad(find(diff(minRad)~=1)+1:size(minRad,1));
        end
        if (sum( isnan(rads(indexs)) ) == length(indexs)  ) %if it has been decided that a segment go falsed labeled a turn but is a straight
            turnPosStart(ii) = 1000;
           
        else
        apexPos = cat( 1,apexPos,minRad );
        apexStart = cat(1,apexStart,minRad(1));
        apexEnd = cat(1,apexEnd,minRad(length(minRad)));
        end
       
end
turnPosStart(turnPosStart==1000)=[];

end

