function [ diffOut ] = unequalDiff( x,y, x0 )
%This function differentates data that is unequally spaced using Lagrange
%polynomal interpolation
%   x - set of 3 x cordinates
%
%   y - set of 3 y cordinates 
%
%   x0 - x location where the differentation should occur

L1 = ( 2.*x0-x(2)-x(3) )./( (x(1)-x(2)).*(x(1)-x(3)));

L2 = ( 2.*x0-x(1)-x(3) )./( (x(2)-x(1)).*(x(2)-x(3)));

L3 = ( 2.*x0-x(1)-x(2) )./( (x(3)-x(1)).*(x(3)-x(2)));

diffOut = y(1).*L1 + y(2).*L2 + y(3).*L3;

if (isnan(diffOut) )
    %diffOut = inf;
    diffOut = 1000;
end

if( diffOut == Inf | diffOut == -Inf )
    diffOut = 1000;
end

end

