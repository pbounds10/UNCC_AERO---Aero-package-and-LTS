function [effRadOut] = effRad(rUnloaded,fz,cfz,d0,B,D,E)
%find the effective radius of the tire
%Inputs:
%   rUnloaded - unloaded rad of tire (m)
%   cfz - vertical stiffness (N/m)
%   d0 - deflections at fz0 (m)
%   B - fit coefficient
%   D - fit coefficient
%STD values for B D and E 
%Radial New = B-10, D-.4, E-0.03
%Radial Worn = B-10, D-.2, E-.03

d = fz/cfz;

defNorm = d/d0;

deformation = d0.*(D.*atan(B.*defNorm)+E*defNorm);

effRadOut = rUnloaded - deformation;
end

