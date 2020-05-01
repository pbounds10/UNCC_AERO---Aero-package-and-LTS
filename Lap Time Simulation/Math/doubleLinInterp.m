function [ yOut ] = doubleLinInterp( x1, x2, xInterp, y1, y2, yInterp, q11, q21, q12, q22 )
%performs a double linear interpolation on table data
%   |----|x1 |  xInt |x2 |
%   |y1  |q11|       |q21|
%   |yInt|   |  yOut |   |
%   |y2  |q21|       |q22|

den = (x2-x1).*(y2-y1);

p1 = ( (x2-xInterp).*(y2-yInterp) )./den.*q11;
p2 = ( (xInterp-x1).*(y2-yInterp) )./den.*q21;
p3 = ( (x2-xInterp).*(yInterp-y1) )./den.*q12;
p4 = ( (xInterp-x1).*(yInterp-y1) )./den.*q22;

yOut = p1 + p2 + p3 + p4;

end

