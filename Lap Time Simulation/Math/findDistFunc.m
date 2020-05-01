function [ xOut, kk ] = findDistFunc( xStart,xGuess,fx,maxDist,minDist,increment )
%Determines the x location that a distance away from point xStart
%Inputs  
%   xStart -  the location that you want the distance from
%   xGuess -  the first x location guess
%   fx - the function that you wish to use
%   maxDist - the maximum distance between points
%   minDist - the minimum distance btween the points
%Outputs
%   xOut - the x location that satisifies the parameters

%Logic
%if the distance is greater than the maxDist then the xGuess must be
%reduced to bring it closer to the xStart
%If the distance is less than the specified minDist then xGuess must be
%increased to get it farther away from xStart
%if the increment of the xGuess is too coarse then you may jump back in
%forth over the solution so a method of increment refinement must be used.
%This saves time over a very find increment as the starting location may be
%far way from the correct location

xGuessInt = xGuess;

%The distance formula using input function
%fDist = @(y) sqrt( (xStart-y).^2 + (fx(xStart)-fx(y)).^2);
fDist = distCalc( xStart,fx(xStart),xGuess,fx(xGuess) ); 

%Check to see if the initial guess is close enough
if(fDist < maxDist && fDist > minDist)
    xOut = xGuess;
    kk = 0;
    return
end    

xGuessOld(1) = 0; %this is setup to make sure that we don't get locked jumping to either side of the answer
kk = 2; %iteration counter
signAug = sign(xGuessInt-xStart);
%inc = signAug.*abs(increment); %the increment size

inc = signAug.*(xGuessInt-xStart)./10;

while(1)
    
    if (kk > 5)
        if( xGuessOld(kk-2) == xGuess ) %if the last guess was exactly the same as this one it means we were juemping around the answer
            inc = inc - (signAug.*0.1); %lower the increment so that we can get to a solution
        end
    end
    
    xGuessOld(kk) = xGuess; %save the old answer
  
    yStart = fx(xStart);
    yGuess = fx(xGuess);
    
    fDist(kk) = round(distCalc( xStart,yStart,xGuess,yGuess ),1);
    
    if ( fDist(kk) <= minDist )
        xGuess = xGuess + inc; %increse the x guess so that the distance can increase
    elseif ( fDist(kk) >= maxDist )
        xGuess = xGuess - inc;
    end
    
    
    if(fDist(kk) <= maxDist & fDist(kk) >= minDist) %If the distance falls between the max and min then use this value
        xOut = xGuess;
        disp(xOut);
        break
    end 
    kk = kk + 1;
    
    if kk > 3000
       fDist = [];
       xGuessOld = [];
       kk = 1;
       xGuess = xGuessInt;
       inc = inc./10;
    end
end    

end

