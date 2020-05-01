function [xNew] = methodBisection(fx,xPos,xNeg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xNew = 1;
while(1)

         xOld = xNew;

         xNew = (xPos+xNeg)./2;

         if (fx(xNew) < 0 )
             xNeg = xNew;
         else
             xPos = xNew;
         end

         if (abs(xNew-xOld) < eps)
             break
         end


    
end

end

