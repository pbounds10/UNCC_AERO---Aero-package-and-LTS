function [ angleOut ] = vectorAngle( vec1, vec2 )
%calculates the angle between two vectors
%Inputs:
%   vec1 - vector 1
%   vec2 - vector 2
%Outputs
%   angleOut - angle betweent the two vectors in degrees

num = vec1(1,1).*vec2(1,1)+vec1(1,2).*vec2(1,2);

den = sqrt( vec1(1,1).^2+vec1(1,2).^2 ).*sqrt( vec2(1,1).^2+vec2(1,2).^2 );

angleOut = acosd(num./den);


end
