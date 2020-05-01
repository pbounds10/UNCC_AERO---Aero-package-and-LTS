function [ diffs ] = smoothData( pos, maxDist )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for ii = 1:length(pos)
    
    if( ii == 1)
    xDiff = pos(ii:ii+2,1);
    yDiff = pos(ii:ii+2,2);
    diffs(ii,1) = unequalDiff(xDiff,yDiff,xDiff(1));
    diffs(ii,2) = unequalDiff(yDiff,xDiff,yDiff(1));
    elseif(ii == length(pos))
    xDiff = pos(ii-2:ii,1);
    yDiff = pos(ii-2:ii,2);
    diffs(ii,1) = unequalDiff(xDiff,yDiff,xDiff(3));
    diffs(ii,2) = unequalDiff(yDiff,xDiff,yDiff(3));
    else
    xDiff = pos(ii-1:ii+1,1);
    yDiff = pos(ii-1:ii+1,2);
    diffs(ii,1) = unequalDiff(xDiff,yDiff,xDiff(2));
    diffs(ii,2) = unequalDiff(yDiff,xDiff,yDiff(2));
    end
    
end

%diffs(abs(diffs) > 20) = 0;

%scatter( pos(:,1), pos(:,2), 20, diffs(ii,2), 'filled' )

end

