function [ vector ] = makeVector( pos1,pos2 )
%makes a vector out of two points
%Inputs
%   pos1 - x and y location of first point
%   pos2 - x and y location of second point
%Outputs
%   vector - 2d vector of x and y


vector(1,1) = pos2(1,1) - pos1(1,1);
vector(1,2) = pos2(1,2) - pos1(1,2);


end

