function [ yOut ] = linInterp( x1, y1, x2, y2, xInterp)
%linearly interpolate between two points
%   Detailed explanation goes here

yOut = y1 + (xInterp - x1).*(y2-y1)./(x2-x1);


end

