function [ pos ] = agumentTrackPos( pos )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for ii = 1:length(pos)
    pos(ii,1) = pos(ii,1) + (-.05 + (.05+.05)*rand(1));
end

end

