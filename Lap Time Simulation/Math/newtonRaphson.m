function [xn1] = newtonRaphson(ft,ftPrime,x0,tol)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

xn1 = x0;

while(1)
    
    xn = xn1;
    
    xn1 = xn - ft(xn)./ftPrime(xn);
    
    err = abs(xn-xn1);
    
    if (err < tol)
        break
    end


end

