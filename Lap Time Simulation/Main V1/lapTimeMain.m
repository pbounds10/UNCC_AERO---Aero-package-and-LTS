clear
clc
%Transmission
gearRatio = [3.285, 1.96, 1.322, 1.028, .820];
finalDrive = 4.176;
rpm = [3500, 4000, 4500, 5000, 5500, 6000];
torque = [111, 113, 113, 108, 102, 92]; %ft-lbf
torque = torque.*1.3558179483314; %in N-m

%BODY DIMENSIONS & WEIGHTS
wheelBase = 2.400; %kg (94.5in)
centerGravity = .762; %m cg hieght(30inches) ...1986 mr2 v153a and v153b: .314 .284 (sae1999-01-1336)
b = 1.3; %m (cg distance from FRONT axle)...1986 mr2 v153a and v153b:1.314 1.284

wfs = 466.293; %kg (1028 lb) 
wrs = 643.6476; %kg (1419 lb) 
mass = wfs + wrs; %kg

%ATMOSPHERE
temperature = 20+273.15; %Kelvin = Degrees C + 273.15
pressure = 101000; %Pa
R = 287.05; %J/kgK
rho = pressure/(R*temperature); %kg/m^3
g = 9.81; %m/s^2

%AERODYNAMICS
cd = .424; % (Drag / X Forces)
cdf = .234; %downforce not lift (Negative is making lift) (Downforce / Z Forces)
area = 1.681; %m^2 (Frontal Area)
kDrag = .5*rho*cd*area;
%aeroCD = %importAero('aeroCos.xlsx','cd','A1:G6');
%aeroCDF = %importAero('aeroCos.xlsx','cdf','A1:G6');
%aeroFront = importAero('aeroCos.xlsx','percentFront','A1:G6');

%WHEELS & TIRES
tireCo = 1; %soon to be removed
tireCos = importTires2("C:\Users\Patrick\Desktop\College Spring 2020\SD 2\Lap Time Sim\Lap Time 3-19-2020\Lap Time 3-19-2020\Lap Time Simulation\Tires\tireCoefficients.xlsx", "Sheet1", [2, 4]);

wheelRadius = 0.31115; %meters
wheelRadiusFront = 0.3149; %meters
wheelRadiusRear = 0.3259; %meters
camberFront = 3;
camberRear = 3;

%BRAKING
    %BIAS
    %used https://www.tceperformanceproducts.com/bias-calculator/ and brake
    %specs to determine brake gain, (rotor torque/psi) per front and rear.
    bgf = 5.825792501; %front brake gain (in-lb/psi)
    bgf = bgf*703.07/86.79616597966708; %(kg-m/(kg/m^2))
    bgr = 4.131717865; %rear brake gain (in-lb/psi)
    bgr = bgr*703.07/86.79616597966708; %(kg-m/(kg/m^2))
    bbias = .567; % percentage of front bias 
    %PROPORTIONING
    %proportioning valve 426.69/.6 (psi/%) = 300000 kg/m^2
    pCritical = 300000; %kg/m^2
    proportion = 0.6; 
    
    %BRAKING PRESSURES
%         pApplied = [0:100:1200]; %Braking Pressure Applied (psi)
%         pApplied = 1000; %psi (CHOOSE MAXIMUM BRAKING PRESSURE)
%         pApplied = pApplied*703.07; %kg/m^2
%         pFront = pApplied;
%         if (pApplied < pCritical)
%             pRear = pApplied;
%         else
%             pRear = pCritical + proportion(pApplied - pCritical);   
%         end

        
    
    accel = 8.5; %braking acceleration m/s^2

%import data REMEBER TO CHANGE THE IMPORT VALUES
%speedAct = importSpeed('winstonProcessed.xlsx','Sheet1',1,952);
[x, y, turn ] = importfile('winstonProcessed.xlsx','Sheet1',1,952);
%[x, y, turn ] = importfile('danvilleProcessedData.xlsx','Sheet1',1,1005);
%[x, y, turn ] =importfile('basic oval track calculations.xlsx','Sheet2',30,281);
%[x, y, turn ] =importfile('ovalWithMorePoints.xlsx','Sheet1',2,662);
x = str2double(x);

%TRACK INITIALIZATION

    %put data into a 2d array
    pos = [x ,y]; % lat & long (m)

    %initialize radius and distance
    %R = zeros(1); %renamed to rads
    d = zeros(1); %distance between all points

    %give a starting velocity (code doesnt like to start from 0)
    vIntial = 12; %m/s
    gear = 1; %the starting gear
    time = 0; %the total time for navigate the track


%%%%%%%%%%%%%%START LAP TIME SIMULATION%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%CHARACTERIZE turns: GET Radius, Start, and End positions%%%%%%%%

    [numTurns, turnPosStart, turnPosEnd, rads] = defineCorners(pos,turn);
    %Start and Ends are determined manually as of 1/16

%this makes a quick initial pass to roughly calculate the velocities but it
%does not give a smooth function aka it doesn't brake... they gon crash

turnTrue = 0;
kk = 1;
v = zeros(0);
v(1) = vIntial;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INITIAL PASS FOR VELOCITIES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ii = 2:length(pos)

             if (turn(ii-1) == 0 && turn(ii) == 1) %check to see if we are in a turn
                turnTrue = 1;
             end

             if (turnTrue == 1)
                 if (ii ~= length(pos))
                    if (ii == 2) %if we happen to start in a curve
                        %disp("in turn " + num2str(ii))
                        v(ii) = curveVelo(mass, rads(kk), area, cd, tireCo, cdf, vIntial);
                    else
                        %disp("in turn " + num2str(ii))
                        v(ii) = curveVelo(mass, rads(kk), area, cd, tireCo, cdf, v(ii-1));
                    end
                    %disp("in curve")
                     if(turn(ii) == 1 && turn(ii+1) == 0) %Check to see if we should leave the turn
                       % disp("off turn " + num2str(ii))
                        turnTrue = 0; % no longer in a turn
                        kk = kk + 1; %proceed to next turn
                     end

                 elseif(ii == length(pos)) %check for end of track
                    
                    [tireMax, tireStiffRear] = magicTire(wRear,1,camberRear,tireCos,3);
                     
                    v(ii) = curveVelo(mass, rads(kk), area, cd, tireCo, cdf, vIntial);
                    %disp("in curve")
                    if(turn(ii) == 1 && turn(1) == 0) %11
                        turnTrue = 0; 
                        kk = kk + 1;
                     end
                 end

             end
                %Find the distance between points
                if (ii < length(pos)) %if we haven't hit the last point then keep on going
                    d(ii) = pdist( [pos(ii,:); pos(ii+1,:)] ,'euclidean');
                else % if we have hit the last point then take the last point and the first point
                    d(ii) = pdist( [pos(ii-1,:); pos(1,:)] ,'euclidean');
                end   

            if (turnTrue == 0) %if there is no turn use the straights calculation 

                if (ii == 2) %if its the first iteration use and intial velocity
                    %disp("in straight")
                    %v(ii) = strLineVelo(d(ii), powEff, vIntial, area, cd, mass);
                    thrust = thrustForce(vIntial, gearRatio, finalDrive, torque, wheelRadius,rpm); %thrust force output from engine in N
%                     wRear = wrs./2.*g./1000; %normal force on one rear tire in kN
%                     
%                     
%                     [tireMax, tireStiffRear] = magicTire(wRear,5,camberRear,tireCos,3);
%                     tireMax = tireMax.*2;
%                     if( thrust > tireMax )
%                         disp('too much thrust from engine')
%                         thrust = tireMax;
%                     end
                    
                    v(ii) = strLineVelo2(d(ii), thrust,vIntial, area, cd, mass);
                else
                    thrust = thrustForce(v(ii-1), gearRatio, finalDrive, torque, wheelRadius,rpm);
                    %[tireMax, tireStiffRear] = magicTire(wRear,5,camberRear,tireCos,3);
%                     tireMax = tireMax.*2; 
%                     if( thrust > tireMax )
%                         disp('too much thrust from engine')
%                         thrust
%                         tireMax
%                         thrust = tireMax;
%                     end
                    v(ii) = strLineVelo2(d(ii), thrust,v(ii-1), area, cd, mass);
                end
            end


        end

        v(length(v)) = v(length(v)-1);
        d(length(v)) = d(length(v)-1);

%%%%%%%%%%%%%%%%%%%%%%%END INITIAL PASS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vp1 = v'; %trouble shooting variable

%%%%%%%%%%%%%%%%%%%DETERMINE ACCELERATION OR BRAKING%%%%%%%%%
%now we will find the max speeds for the turns

vTurn = zeros(0); %this is where we will store our max turn velocity
kk = 1;

turnAorD = zeros(0); %1 means when need to accelerate in the turn and 0 means we need to break

%find the max turn velocity by find the velocity at the first point in the
%turn
for ii = 2:length(pos)
    if (turn(ii-1) == 0 && turn(ii) == 1)  %check to see if we are in a turn 11
        turnTrue = 1;
    end
    
    if (turnTrue == 1)
        vTurn(kk) = v(ii); %when we come to a turn save the stored velocity
       
        if (vTurn(kk) > v(ii-1))
            turnAorD(kk) = 1; %accelerating
        else
            turnAorD(kk) = 0; %braking
        end
        kk = kk + 1; 
        turnTrue = 0;
    end       
end
%%%%%%%%%%%%%%%%%%%END DETERMINE ACCELERATION OR BRAKING%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%MAKE VELOCITIES CONTINOUS WITH BRAKING OR ACCELERATION%%%%%%%%%%%%%%%

brakeSpeed = zeros(0);%this is used if we need to brake
accelSpeed = zeros(0);%this is used if we need to speed up in the turn

%This will create a continous velocite around the track
for bb = 1:length(turnPosStart) %This walks through the curves
  
    if (turnAorD(bb) == 1) %this means that we need to accelerate in the turn
        for ii = turnPosStart(bb):1:turnPosEnd(bb)
            v(ii) = curveAccelVelo( mass, rads(bb), area,cd,tireCo,cdf,v(ii-1), d(ii) );
        end     
    end
    
    
    
    if (turnAorD(bb) == 0) %if we need to break do this section of code
    brakeSpeed = zeros(0);
    
    if (bb == 1) %check out if this check in neeeded
        
        for ii = turnPosStart(bb):-1:1 %if we have hit the end of the track do this
           if (ii == turnPosStart(bb))  %if we are at the start position use the max Turn velocity
                brakeSpeed(ii) = constABraking(accel, v(ii), d(ii));
                %*** HERE***
           else
                brakeSpeed(ii) = constABraking(accel, brakeSpeed(ii+1), d(ii)); %it's plus 1 because we are going backwards  
           end

        end
        
        [ind, err1] = straightBrake(v,brakeSpeed, turnPosStart(bb),1 );
        for cc = turnPosStart(bb):-1:ind
           v(cc) = brakeSpeed(cc);
              
        end
    else
       for ii = turnPosStart(bb):-1:1 %turnPosEnd(bb - 1) %if we havent hit the end of the track do this
           if (ii == turnPosStart(bb)) %if we are at the start position use the max Turn velocity
                brakeSpeed(ii) = constABraking(accel, v(ii), d(ii));
           else
                brakeSpeed(ii) = constABraking(accel, brakeSpeed(ii+1), d(ii)); %it's plus 1 because we are going backwards  
           end
       end 


       [ind, err1] = straightBrake(v, brakeSpeed, turnPosStart(bb), 1); %, d, accel, bb );

       for cc = turnPosStart(bb):-1:ind
           v(cc) = brakeSpeed(cc);
           
       end
    end 
    
    end %end of braking section
end

v(length(v)) = v(length(v)-1);
d(length(v)) = d(length(v)-1);

%%%%%%%%%%%%END  MAKE VELOTIES CONTINOUS WITH BRAKING OR ACCELERATION%%%%%%%%%%%%%%%

%%%%%%%%%DETERMINE SHIFTING TIME%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0;

gearPos = [];
bb = 1;
gearShiftNum = 0;

for ii = v
    gearPos(bb) = gear;
    [dt, gear] = gearShift( gear, ii, gearRatio, finalDrive , wheelRadius, rpm);
    bb = bb + 1;
    %time = time + dt;
    
    if (dt ~= 0 )
        gearShiftNum = gearShiftNum + 1;
    end
end

disp("Shifting added " + num2str(time) + " seconds")

gearPos = gearPos';


%%%%%%%%%%%%%%%%DETERMINE TIME FROM VELOCITY FIELD%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0;

for ii = 1:1:length(pos)
    if ( ii == 1)
        vAvg = v(ii);
        dt = d(ii)/vAvg;
    else
        vAvg = ( v(ii)+v(ii-1) )/2;
        dt = d(ii)/vAvg;
    end
    timeAccum(ii) = time;
    time = time + dt;
    
    
end
 
disp("The final time is: " + num2str(time))
 
speedMPH = v.*2.237;

speedMPH = speedMPH';

% deltaV = speedAct - speedMPH; 

d = d';
v = v';
pointsize = 30;

figure(1)
hold
plot(1:size(pos),speedMPH, 'Linewidth',4)
% plot(1:size(pos),speedAct, 'Linewidth',4)
xlabel('Track Position')
ylabel('Speed (MPH)')
legend('Predicted Speed', 'Actual Speed')


% scatter(x, y, pointsize, v*2.237, 'filled');
% title("Velocity (Mph)")
% xlabel("meters")
% ylabel("meters")
% figure(2)
% scatter(x, y, pointsize, turn, 'filled');
% title("Danville Autocross Track")
% xlabel("meters")
% ylabel("meters")
% figure(3)
% scatter(x, y, pointsize, deltaV, 'filled');
% title("Delta Velocity (Mph)")
% xlabel("meters")
% ylabel("meters")


%figure(2)
%scatter(x, y, pointsize, d, 'filled');
