function [ distOut ] = distCalc( x1,y1,x2,y2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

distOut = sqrt( (x1-x2).^2 + (y1-y2).^2 );

end

