function [xn1] = NewtonRapOpt(fp,fpp,x1, tol)
%Finds the maxima of a functin using the 1st and 2nd derivative
%   Detailed explanation goes here

xn = x1; %initialize the values
k = 1; %iteration counter

while(1)
    
    xn1 = xn - (fp(xn)./fpp(xn));
    
    dx = xn1 - xn;
    
    if (abs(dx) < tol)
        break
    end
    
    k = k + 1;
    xn = xn1;
    
end
    


end

