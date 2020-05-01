function [ indexLocations ] = findInflection( inputArray )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

indexLocations = [];

checkLength = 4;

for ii = 1:length( inputArray )
    
    aa = ii - 1;
    bb = ii + 1;
    
    if (ii <= checkLength)
        if ( abs( inputArray(ii) ) > abs( inputArray(1:aa) ) & abs( inputArray(ii) ) > abs( inputArray(bb:ii+checkLength) ) )
            indexLocations = cat(1,indexLocations,ii);
        end
    elseif( ii > length(inputArray) - checkLength )
        if ( abs( inputArray(ii) ) > abs( inputArray(ii-checkLength:aa) ) & abs( inputArray(ii) ) > abs( inputArray(bb:length(inputArray)) ) )
            indexLocations = cat(1,indexLocations,ii);
        end
    else
        if ( abs( inputArray(ii) ) > abs( inputArray(ii-checkLength:aa) ) & abs( inputArray(ii) ) > abs( inputArray(bb:ii+checkLength) ) )
            indexLocations = cat(1,indexLocations,ii);
        end
    end


end

end

