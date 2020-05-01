function [ R ] = radFinder(p1,p2,p3)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

aa = p1(1,:);
bb = p2(1,:);
cc = p3(1,:);

cDis = ((bb(1)-aa(1))^2+(bb(2)-aa(2))^2)^.5;
bDis = ((cc(1)-aa(1))^2+(cc(2)-aa(2))^2)^.5;
aDis = ((cc(1)-bb(1))^2+(cc(2)-bb(2))^2)^.5;

inter = (aDis^2+bDis^2+cDis^2)/(2*bDis*cDis) - 2
angle_A = acosd( inter );

angle_P = 180 - angle_A;

R = aDis/(2*sind(angle_P));

end

