function [ turnOutput ] = findCurveStart( inflections, slopes )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

checkLength = 4;

turnOutput = [];

indexLocation = [];


for ii = 1:length(inflections)
    err1 = false;
    

    for aa = inflections(ii)-1:-1:inflections(ii)-20
        
        bb = aa-1;
        
        if ( aa - checkLength <= inflections(ii)-20 )
            inflections(ii) = 1000;
            err1 = true;
            continue
        end

        
%         aa-checkLength:aa-1
%         
%         abs( slopes(aa) )- abs( slopes(aa-checkLength:aa-1)) > 1 
%         
%         abs( slopes(aa-checkLength)) - abs( slopes(aa-checkLength-1) ) < .5
%         
        slopes(aa-checkLength:aa)
        aa
        ( abs( slopes(aa) ) - abs( slopes(bb) ) > 1 )
        ( abs( slopes(bb) ) - abs( slopes(bb-checkLength:bb-1) ) < 1.5 )
        
        ( abs( slopes(aa) ) - abs( slopes(bb) ) > 1 ) & ( abs( slopes(bb) ) - abs( slopes(bb-checkLength:bb-1) ) < 1.5 )
        
        if ( ( abs( slopes(aa) ) - abs( slopes(bb) ) > 1 ) & ( abs( slopes(bb) ) - abs( slopes(bb-checkLength:bb-1) ) < 1.5 ) )
            indexLocation = aa;
            disp(['the turn starts at ', num2str(indexLocation)])
            break
        end
        
        
        
    end
    
    if (err1 == false)
        turnOutput = cat(2,turnOutput, indexLocation:inflections(ii) );
    end
    
    
    
    
end
        
 inflections   

end

