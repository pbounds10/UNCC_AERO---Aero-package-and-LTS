clear
clc
%[Latitude, Longitude] = importRawTrack("M:\projects\Dhillon-00AEROUNCCSD1\lap time sim 3\Log-20190825-104416 Danville 6 - 0.40.710 - Interpolated.csv", [4934, Inf]);
[Latitude, Longitude,speed] = importRawTrackwSpeed('Log-20191020-134416 winston 8.19.20 - 0.45.949 - Interpolated.csv',11050, 12001);
%[Time,UTCTime,Latitude,Longitude,Altitudem,Altitudeft,SpeedMPH] = importfile1('Log-20200116-212751 Data Log - Interpolated.csv',5, 73);


Latitude = rmmissing(Latitude);
Longitude = rmmissing(Longitude);


kmLat = Latitude.*(10000/90)*1000; %converts lat degrees to meters
kmLong = Longitude.*(10000/90)*1000; %converts long degrees to meters

for ii = 1:length(kmLat)
   
     y(ii) = round( kmLat(ii) - kmLat(1), 2);
     x(ii) = round( kmLong(ii) - kmLong(1), 2);   
    
end
y = y';
x = x';
pos = horzcat(x,y);

% danville posFind = [0,0; 46.04,-40.46; -109.02, -143.83; -166.11,-169.79; -232.28,-205.61; -199.4,-225.69; 26.91, -100.87; 74.19, -73.18; 78.48,-72.78; 110.03,-53.32; 204.18,-11.22; 173.89,10.46;];
%Winston
%posFind = [102.12,4.2; 107.53,10.5; 108.14,14.51; 96.86, 20.24; 80.68,6.24; 78.44,-4.58;  75.92, -45.52; 76.13, -76.10; 86.56, -125.17; 111.24,-143.26; 160.21,-114.32; 158.86,-78.39; 158.22, 29.44; 159.56,75.37; 154.57,83.31; 153.96, 101.09; 153.57,142.74; 140.98 ,167.17; 106.81,169.96; 83.01,150.22; 80.28, 141.57; 66.9,126.03; 60.51,122.29; 48.63,90.16 ]; 

clickPos = [];

while(1)

scatter(x,y)
[xTemp, yTemp] = ginput;
posTemp = horzcat(xTemp,yTemp);

clickPos = cat(1,clickPos,posTemp);


    if( mod(size(clickPos,1),2) == 0 )
        break
    else
        warning('Must select an even number of points. Previous points have been saved. Simplely add an end point to the last selection')
    end

end

%posFind = [76.8,-7.6; 106.3,8.5; 104.8,19.1; 84.9,14.4; 78.8,-1.2; 80.8,-21.3; 78.1,-40.4; 71.6,-60.4; 80.9,-86.2; 83.6,-106.8; 92.6,-138.2; 122.8,-142.7; 166,-98.3; 158.3,-77; 154.4,-59.5; 156.9,-38.2; 154.9,-16.2; 157.1,-0.6; 157.7,10.4; 157,24.3; 163.8,38.1; 167.5,53.6; 152.5,91.3; 155.1,106.5; 152.3,134.4; 154,147.8; 147.7,162.6; 122,171.9; 86.4,155.7; 78.5,136; 54.8,116.9; 48.2,101.7; 48.9,78.4; 45.3,65.2; 31.4,40.7];

curve = zeros(size(pos,1),1);

lenPos = size(pos,1);

for ii = 1:2:size(clickPos,1)-1
    %index = find([posFind(ii,1),posFind(ii,2)] == pos(:) )
    %index = find(posFind(ii,1) == pos(:,1) & posFind(ii,2)==pos(:,2));
    for aa = 1:size(pos,1)  
        distStart(aa) = distCalc(pos(aa,1),pos(aa,2),clickPos(ii,1),clickPos(ii,2));
        distEnd(aa) = distCalc(pos(aa,1),pos(aa,2),clickPos(ii+1,1),clickPos(ii+1,2));
    end
    
    indexStart = find( distStart == min(distStart) );
    indexEnd = find( distEnd == min(distEnd) );
    
    curve(indexStart:indexEnd) = 1

%     for bb = 1:index
%         if(bb < lenPos)
%             if ((pos(bb,1) == posFind(ii,1) ) && (pos(bb,2) == posFind(ii,2) ) )
%                 %pos(index,1)
%                 %pos(index,2)
%                 curve(ii) = bb;
%             end
%         end
% 
%     end
end

% for ii = 1:size(posFind,1)
%     %index = find([posFind(ii,1),posFind(ii,2)] == pos(:) )
%     index = find(posFind(ii,1) == pos(:,1) & posFind(ii,2)==pos(:,2));
% 
% 
%     for bb = 1:index
%         if(bb < lenPos)
%             if ((pos(bb,1) == posFind(ii,1) ) && (pos(bb,2) == posFind(ii,2) ) )
%                 %pos(index,1)
%                 %pos(index,2)
%                 curve(ii) = bb;
%             end
%         end
% 
%     end
% end


% turnTrue = 0;
% kk = 1;
% 
% for ii = 1:lenPos
%     if (turnTrue == 0)
%         %disp("true")
%         if(ii == curve(kk) )
%             %disp('turn on')
%             turnTrue = 1;
%             kk = kk + 1;
%             if (kk > length(curve))
%                 kk = length(curve);
%             end
%         end
%     elseif( turnTrue == 1)
%         if(ii == curve(kk) )
%             %disp('turn off')
%             turnTrue = 0;
%             kk = kk + 1;
%             if (kk > length(curve))
%                 kk = length(curve);
%             end
%         end 
%     end
%     
%     if(turnTrue == 1)
%         turn(ii) = 1;
%     else
%         turn(ii) = 0;
%     end
% 
%     
%     
% end
% turn = zeros(length(pos));
% turn = turn';
% 
% dataExport = horzcat(x,y,turn);
% scatter(x,y, 10)
% scatter(x,y, 10, turn)
% 
% %curve = [0; 89; ];
% %index = find([45.1444,-42.6333] == pos(:) )