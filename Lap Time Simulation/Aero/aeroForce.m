function [ frontDF, rearDF, drag ] = aeroForce(yaw, pitch, vX, vY, aeroCD, aeroCDF, aeroFront, area, density, symetric )
%calculate the aerodynamic force on the vehicle
%Inputs
%   yaw - degrees of yaw (only positive if vehicle symetric)
%   pitch - degrees of pitch
%   vX - velocity in x direction
%   vY - velocity in y direction
%   aeroCD - table of cd values
%   aeroCDF - table of cdf values
%   aeroFront - table of %Front values
%   area - frontal area in m^2
%   density - density of air kg/m^3
%   symetric - is the vehicle symetric 1=yes,0=no
%Outputs:
%   frontDF - front downforce in Newtons
%   rearDF - read downforce in Newtons
%   drag - drag force in Newtons
%Logic:
%All aero values will be passed in through a pitch/yaw table. each value
%will then be interpolated through linear double interpolation. the cd,cdf,
%and %Front values will then be used to calculate aerodynamic forces


%if the vehicle is symetric then there is no difference in pos/neg yaw
if ( symetric == 1 )
    yaw = abs(yaw);
end

%find where the yaw fits into the matrix
x1 = find( aeroCD(1,:) <= yaw,1,'last');
x2 = find( aeroCD(1,:) > yaw,1);
    if ( isempty(x2) )%if the yaw is the last value or larger than the matrix values
        x2 = x1;
        x1 = x2 - 1;
    end

%find where the pitch fits into the matrix
y1 = find( aeroCD(:,1) <= pitch,1,'last');
y2 = find( aeroCD(:,1) > pitch,1);
    if ( isempty(y2) )%if the pitch is the last value or larger than the matrix values
        y2 = y1;
        y1 = y2 - 1;
    end

%perform the interpolation for the aero coefficients
cd  = doubleLinInterp(aeroCD(1,x1),aeroCD(1,x2),yaw,aeroCD(y1,1),aeroCD(y2,1),pitch,aeroCD(y1,x1),aeroCD(y1,x2),aeroCD(y2,x1),aeroCD(y2,x2));
cdf = doubleLinInterp(aeroCDF(1,x1),aeroCDF(1,x2),yaw,aeroCDF(y1,1),aeroCDF(y2,1),pitch,aeroCDF(y1,x1),aeroCDF(y1,x2),aeroCDF(y2,x1),aeroCDF(y2,x2));
aeroFront = doubleLinInterp(aeroFront(1,x1),aeroFront(1,x2),yaw,aeroFront(y1,1),aeroFront(y2,1),pitch,aeroFront(y1,x1),aeroFront(y1,x2),aeroFront(y2,x1),aeroFront(y2,x2));

%calcuate the resolved velocity
v = sqrt(vX.^2+vY.^2);

%calcuate the downforce for the whole car
downForce = 0.5.*density.*area.*v.^2.*cdf;

%split the downforce into front and rear
frontDF = downForce.*aeroFront;
rearDF = downForce-frontDF;
%calculate the drag for the whole car
drag = 0.5.*density.*area.*v.^2.*cd;

end

