function [ deltaV ] = deltaSpeed(preVelo,actVelo )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

predV = preVelo.*2.237;

deltaV = actVelo - predV;


end

