clear
clc
%[Latitude, Longitude] = importRawTrack("M:\projects\Dhillon-00AEROUNCCSD1\lap time sim 3\Log-20190825-104416 Danville 6 - 0.40.710 - Interpolated.csv", [4934, Inf]);
[Latitude, Longitude,speed] = importRawTrackwSpeed('Log-20191020-134416 winston 8.19.20 - 0.45.949 - Interpolated.csv',11050, 12001);
%[Time,UTCTime,Latitude,Longitude,Altitudem,Altitudeft,SpeedMPH] = importfile1('Log-20200116-212751 Data Log - Interpolated.csv',5, 73);

Latitude = rmmissing(Latitude);
Longitude = rmmissing(Longitude);

kmLat = Latitude.*(10000/90)*1000; %converts lat degrees to meters
kmLong = Longitude.*(10000/90)*1000; %converts long degrees to meters

%put the data into a 2d position matrix
for ii = 1:length(kmLat)  
     pos(ii,2) = round( kmLat(ii) - kmLat(1), 2); 
     pos(ii,1) = round( kmLong(ii) - kmLong(1), 2);   
    
end


clickPos = [];

%select the data points to turn into segments. Must be an even amount of
%points
while(1)

scatter(pos(:,1),pos(:,2)) %show the graph of the track
title("Select and even number of points. Press ''Enter'' when done selecting points")

[xTemp, yTemp] = ginput; %select the data points you want
posTemp = horzcat(xTemp,yTemp);

clickPos = cat(1,clickPos,posTemp); %add the points to the matrix

    if( mod(size(clickPos,1),2) == 0 ) %if the user has selected an even number of points continue
        break
    else
        warning('Must select an even number of points. Previous points have been saved. Simplely add an end point to the last selection')
    end

end


turn = zeros(size(pos,1),1);

lenPos = size(pos,1); 

indCount = 1; %store all the indexs for later use of saving

for ii = 1:2:size(clickPos,1)-1

    for aa = 1:size(pos,1)  
        distStart(aa) = distCalc(pos(aa,1),pos(aa,2),clickPos(ii,1),clickPos(ii,2)); %find the distance between the odd point and the rest of the points
        distEnd(aa) = distCalc(pos(aa,1),pos(aa,2),clickPos(ii+1,1),clickPos(ii+1,2)); %find the distance between the even point and the rest of the points
    end
    
    indexStart(indCount) = find( distStart == min(distStart) ); %find the index location of the start of the segment
    indexEnd(indCount) = find( distEnd == min(distEnd) ); %find the index location of the end of the segment
    
%     if (indexStart(indCount) < indexStart(1:indCount))
%         
%     end       
%         
%     if (indexEnd(indCount) < indexEnd(1:indCount))
%         
%     end
    
    turn(indexStart(indCount):indexEnd(indCount)) = 1;
    
    indCount = indCount + 1;

end

dataExport = horzcat(pos(:,1),pos(:,2),turn);
% scatter(x,y, 10)
scatter(pos(:,1),pos(:,2), 10, turn)
