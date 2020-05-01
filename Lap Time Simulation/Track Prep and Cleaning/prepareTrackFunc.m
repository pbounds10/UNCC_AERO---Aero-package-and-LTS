function [pos, segment, index] = prepareTrackFunc(dataFileIn,row1,row2, dataFileOut, segmentIndex, convertToM)
%This function converts lat and long data to (x,y) meters and  segments track
%INPUTS:
%   row1 - start of data row
%   row2 - end of data row
%   dataFileIn - file name and directory of the data
%   dataFileOut - where to write the data file is desired .xlsx. set as nan
%   to not write to file
%   segmentIndexOld - put 0 if wanting to create new segments.
%   convertToM -  0 = lat and long, 1 = meter input
%OUTPUTS:
%   pos - x y position in meters
%   turn - segmented or not (1 or 0)
%   index - index locations of segements

dataIn = xlsread(dataFileIn,strcat("A",num2str(row1),":","B",num2str(row2)));
if convertToM == 0
    Latitude = dataIn(:,2)
    Longitude = dataIn(:,1)


    Latitude = rmmissing(Latitude);
    Longitude = rmmissing(Longitude);

    kmLat = Latitude.*(10000/90)*1000; %converts lat degrees to meters
    kmLong = Longitude.*(10000/90)*1000; %converts long degrees to meters

    %put the data into a 2d position matrix
    for ii = 1:length(kmLat)  
         pos(ii,2) = round( kmLat(ii) - kmLat(1), 2); 
         pos(ii,1) = round( kmLong(ii) - kmLong(1), 2);    
    end
else
    pos(:,2) = dataIn(:,2); 
    pos(:,1) = dataIn(:,1);
end


clickPos = [];
segment = zeros(size(pos,1),1); %intialize segment with zeros

if ( segmentIndex ~= 0 ) %if previous index data has been given then use that and just add to it
    indCount = length(segmentIndex) + 1;
    index = segmentIndex; 
    for ii = 1:length(index)-1 %use the old index data to segment the track before adding new
        segment(index(ii):index(ii+1)) = 1;
    end       
else
    indCount = 1; %store all the indexs for later use of saving
    index = [];
end


%select the data points to turn into segments. Must be an even amount of
%points

while(1) %segment the track 2 clicks at a time

    segmentOld = segment;
    indexOld = index;
    
% while(1)

scatter(pos(:,1),pos(:,2), 10, segment) %show the graph of the track
title("Select 2 points")

[xTemp, yTemp] = ginput(2); %select the data points you want
clickPos = horzcat(xTemp,yTemp);

%clickPos = cat(1,clickPos,posTemp); %add the points to the matrix

%     if( mod(size(clickPos,1),2) == 0 ) %if the user has selected an even number of points continue
%         break
%     else
%         warning('Must select an even number of points. Previous points have been saved. Simplely add an end point to the last selection')
%     end

% end

for ii = 1:2:size(clickPos,1)-1

    for aa = 1:size(pos,1)  
        distStart(aa) = distCalc(pos(aa,1),pos(aa,2),clickPos(ii,1),clickPos(ii,2)); %find the distance between the odd point and the rest of the points
        distEnd(aa) = distCalc(pos(aa,1),pos(aa,2),clickPos(ii+1,1),clickPos(ii+1,2)); %find the distance between the even point and the rest of the points
    end
    

    index(indCount) = find( distStart == min(distStart) ); %find the index location of the start of the segment
    index(indCount+1) = find( distEnd == min(distEnd) ); %find the index location of the end of the segment
    
%     if (abs( index(indCount)-index(indCount-1)) < 3)
%         index(indCount) = index(indCount)+2;
%         
    
    segment(index(indCount):index(indCount+1)) = 1;
    
    indCount = indCount + 2;

end


scatter(pos(:,1),pos(:,2), 10, segment)

%ask the user if they would like to continue segmenting
userAns = questdlg('would you like to continue?','SegmentTrack', 'Yes', 'No','Remove Last','Yes');

switch userAns
    
    case 'Yes'
        %Continue Loop
    case 'No' %exit the while loop
        break 
    case 'Remove Last' %revert to the last time the user selected data
        index = indexOld;
        segment = segmentOld;
        scatter(pos(:,1),pos(:,2), 10, segment)
end


end

diffs = find( abs(diff(segment)) == 1);

%come back to this and write a rule to make sure that every points had at
%least two consecutive points of 0 or 1
% for ii = 2:length(diffs)
%     if( diffs(ii-1)-diff(ii) < 3)
%         if (diffs(ii-1) > 
%         segment(diffs(ii-1):diffs(ii)+3) = 0;

dataExport = horzcat(pos(:,1),pos(:,2),segment);


if (~isnan( dataFileOut ) )
    fileName = strcat(dataFileOut,".xlsx");
    writematrix(dataExport,fileName,'Sheet',1);
end

end
