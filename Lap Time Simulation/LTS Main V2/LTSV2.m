%LTS V2.0

%This LTS is to incorperate a moving frame of reference and new braking,
%tire, and corner models


%%%%%%%%%% START TEMP %%%%%%%%
addpath(genpath('D:\Desktop\College Spring 2020\SD 2\Lap Time Sim\Lap Time 3-19-2020\Lap Time 3-19-2020\Lap Time Simulation'))
clear
clc
close all
%% import the track data
[x, y, turn ] = importfile('tracks\Cleaned\winstonSegmented.xlsx','Sheet1',1,952);
%[x, y, turn ] = importfile('tracks\Cleaned\danvilleSegmented.xlsx','Sheet1',1,1005);
x = str2double(x);
%% ATMOSPHERE
temperature = 20+273.15; %Kelvin = Degrees C + 273.15
pressure = 101000; %Pa
R = 287.05; %J/kgK
rho = pressure/(R*temperature); %kg/m^3
g = 9.81; %m/s^2

%% Transmission
gearRatio = [3.285, 1.96, 1.322, 1.028, .820];
finalDrive = 4.176;
rpm = [3500, 4000, 4500, 5000, 5500, 6000];
torque = [111, 113, 113, 108, 102, 92]; %ft-lbf
torque = torque.*1.3558179483314; %in N-m

driveWheels = 1;

%% BODY DIMENSIONS & WEIGHTS
wheelBase = 2.400; %kg (94.5in)
centerGravity = .762; %m cg hieght(30inches) ...1986 mr2 v153a and v153b: .314 .284 (sae1999-01-1336)
lenA = 1.3; %m (cg distance from FRONT axle)...1986 mr2 v153a and v153b:1.314 1.284
lenB = wheelBase-lenA;

mfs = 466.293; %kg (1028 lb) 
mrs = 643.6476; %kg (1419 lb) 
mass = mfs + mrs; %kg

wfs = mfs.*g;
wrs = mrs.*g;



%% AERODYNAMICS
area = 1.681; %m^2 (Frontal Area)
aeroCD = importAero('aeroCos.xlsx','cd','A1:G6');
aeroCDF = importAero('aeroCos.xlsx','cdf','A1:G6');
aeroFront = importAero('aeroCos.xlsx','percentFront','A1:G6');

%% WHEELS & TIRES
tireCos = importTires2("Tires\tireCoefficients.xlsx", "Sheet1", [2, 4]);

wheelRadFront = 0.3149; %meters
wheelRadRear = 0.3259; %meters
camberFront = -3; %in degrees
camberRear = -3; %in degrees

%% Spring Rates
springF = 5.*1000; %kg/mm
springR = 8.*1000; %kg/mm

%% BRAKING
    %BIAS
    %specs to determine brake gain, (rotor torque/psi) per front and rear.
    bgf = 5.825792501; %front brake gain (in-lb/psi)
    bgf = bgf*0.1129/6.899476; %(N-m/kpa)
    bgr = 4.131717865; %rear brake gain (in-lb/psi)
    bgr = bgr*0.1129/6.899476; %(N-m/kpa)
    bbias = .567; % percentage of front bias 
    %PROPORTIONING
    %proportioning valve 426.69/.6 (psi/%) = 300000 kg/m^2
    pCritical = 2941.995009; %kg/m^2
    proportion = 0.6; 

    pIn = 5000; %driver input pressure kpa
    pressMax = 8273.709; %kpa
    
%%%%%%%%%%%%%%%% END TEMP %%%%%%%




%% TRACK INITIALIZATION

%put data into a 2d array
pos = [x ,y]; % lat & long (m)

%initialize distance
d = zeros(1); %distance between all points

%sweepVar
alphaSweep = 0:.1:15;
%% intialize velocities and key params
currentPos = 0; %position in track matrix
u = zeros(1,size(pos,1));
v = zeros(1,size(pos,1));
ay = zeros(1,size(pos,1));
ax = zeros(1,size(pos,1));
wF = zeros(1,size(pos,1));
wR = zeros(1,size(pos,1));
alphaF = zeros(1,size(pos,1));
alphaR = zeros(1,size(pos,1));
steerAngle = zeros(1,size(pos,1));

uBack = zeros(1,size(pos,1));
vBack = zeros(1,size(pos,1));
axBack = zeros(1,size(pos,1));
ayBack = zeros(1,size(pos,1));
wFBack = zeros(1,size(pos,1));
wRBack = zeros(1,size(pos,1));

u(1) = 12;
v(1) = 0;
ay(1) = 0;
ax(1) = 0;
alphaF(1) = 0;
alphaR(1) = 0;
fyFront(1) = 0;
fyRear(1) = 0;
wF(1) = wfs; %add downforce
wR(1) = wrs; %add downforce

gear = 1; %the starting gear
time = 0; %the total time for navigate the track
dt = 0;




%% %%%%%%%%%%%%%%%%%%%START LAP TIME SIMULATION%%%%%%%%%%%%%%%%%%%%%%%%%%

%% find distance bewteen points
for ii = 1:size(pos,1)
        %Find the distance between points
    if (ii < length(pos)) %if we haven't hit the last point then keep on going
        d(ii) = pdist( [pos(ii,:); pos(ii+1,:)] ,'euclidean');
    else % if we have hit the last point then take the last point and the first point
        d(ii) = pdist( [pos(ii-1,:); pos(1,:)] ,'euclidean');
    end   
end

for ii = 2:size(pos,1)
    if(ii ~= size(pos,1))
        if (d(ii)>5)
            d(ii) = (d(ii-1)+d(ii+1))./2;
        end
    else
        if (d(ii)>10)
            d(ii) = (d(ii-1)+d(ii-2))./2;
        end
    end
end
        

%% CHARACTERIZE turns: GET Radius, Start, and End positions and Apex Info%%%%%%%%

[turnStart,apexStart,apexEnd,rads,apexInd] = defineCorners(pos,turn);

% for ii = 1:size(apexPos,1)
%     apexInd(ii) = find(apexPos(ii,1) == pos(:,1) & apexPos(ii,2) == pos(:,2));
% end

%% find Apex speeds
for ii = apexInd
    drag = 0;
    [uApex, ayApex, alphaFApex,alphaRApex,steerAngleApex,fyR,fyF] = findLatAccelMax(tireCos,camberFront,camberRear,rads(apexInd),lenA,lenB,wfs,wrs,drag,1);
    u(ii) = uApex;
    v(ii) = 0;
    ay(ii) = ayApex;
    ax(ii) = 0;
    alphaF(ii) = alphaFApex;
    alphaR(ii) = alphaRApex;
    fyFront(ii) = fyF;
    fyRear(ii) = fyR;
    wF(ii) = wfs; %add downforce
    wR(ii) = wrs; %add downforce
end
%% initialize starting point
apexInd = cat(1,1,apexInd);%add the start of the track to the 'apexes"
apexEnd = cat(1,1,apexEnd);
apexStart = cat(1,1,apexStart);


uFinal = u;


%must gen lat forces different because magicTire cannot take arrays

    
%calculate all accels out of apex up to next apex. then go back and
%calculate (back integrate) to where the braking zone meets the
%acceleration zone. trim the data set. maybe try to do this with time
%integration. or just keep with the method I have been using.

%% forward acceleration
%known at apex
%forward acceleration can happen in two different places, in a curve and
%straights
%for staights v = 0, ay = 0; nothing happens in lateral

%for curves:
%rDot = 0; Fyfront*a = fyrear*b vDot = 0(slidding lateral accel)
%m(uDot-vr) = sum(Fx) - v is known
%+mur = sum(Fy)
%assume we are trying to maximize long accel
%lat accel will just to complete the curve. After that all goes to long

for ii = 1:size(apexEnd,1) %do this for every apex
    posInc = 1 + apexEnd(ii);
    %if (ii ~= size(apexStart,1))
   if (ii ~= size(apexStart,1)) %make sure we arent at the end of the track
       endCond = apexStart(ii+1);
   else
       endCond = size(pos,1);
   end%make sure we arent at the end of the track so we dont have an index exception
        while(posInc < endCond)%apexStart(ii+1))%turnStart(ii+1)) %start forward acceleration as long as we haven't reached next apex or end of track
                wF(posInc) = wF(posInc-1);
                wR(posInc) = wR(posInc-1);
            %% if we are in a turn
            if (~isnan(rads(posInc))) 
                %find velocity at new point based on last point data
                u(posInc) = sqrt(u(posInc-1).^2+2.*ax(posInc-1).*abs(pos(posInc,1)-pos(posInc-1,1)));
                v(posInc) = sqrt(v(posInc-1).^2+2.*ay(posInc-1).*abs(pos(posInc,2)-pos(posInc-1,2)));
                
                %curve Radius
                curveRad = rads(posInc);
                %find force balence
                SA = 0;

                ax(posInc) = 1000; %initialize value for use in iterative method
                kappa = 0; %initialize value for use in iterative method
                %% find tire balence
                while(1)
                    axOld = ax(posInc);
                    
                    [~,~,alphaPeakF,~] = magicTire(wF(posInc)/2000,0,0,camberFront,tireCos,1);
                    [~,~,alphaPeakR,~] = magicTire(wR(posInc)/2000,0,0,camberRear,tireCos,1);
                   % find the slip angles of the tires so that traction can
                   % be checked
                    fyFront(posInc) = mass.*(u(posInc).^2./curveRad)./(1+lenA./lenB);
                    fYF_diff = 2.*magicTire(wF(posInc)/2000,alphaSweep,0,camberFront,tireCos,1)-fyFront(posInc);
                    alphaF(posInc) = alphaSweep( find(abs(fYF_diff) == min(abs(fYF_diff))) );
                    if (alphaF(posInc)>alphaPeakF)%if the front tires cant generate the tracion to make the turn
                        alphaF(posInc) = alphaPeakF;
                        fyFront(posInc) = 2.*magicTire(wF(posInc)/2000,alphaF(posInc),0,camberFront,tireCos,1);
                        u(posInc)  = sqrt(fyFront(posInc).*(1+lenA./lenB).*curveRad./mass);
                    end
                        
                    fyRear(posInc) = lenA./lenB.*fyFront(posInc);
                    fYR_diff = 2.*magicTire(wR(posInc)/2000,alphaSweep,kappa,camberRear,tireCos,1)-fyRear(posInc);
                    alphaR(posInc) = alphaSweep( find(abs(fYR_diff) == min(abs(fYR_diff))) );
                   
                   if(driveWheels == 1)  
                        [ax(posInc),wF(posInc),wR(posInc),kappa] = longAccel(u(posInc),v(posInc),gearRatio,finalDrive,torque,wheelRadRear,rpm,wF(posInc),wR(posInc),centerGravity,(lenA+lenB),lenA,alphaR,tireCos,camberRear,driveWheels,springF,springR,0,0,aeroCD,aeroCDF,aeroFront,area,rho,1,curveRad);
                   else
                        [ax(posInc),wF(posInc),wR(posInc),kappa] = longAccel(u(posInc),v(posInc),gearRatio,finalDrive,torque,wheelRadRear,rpm,wF(posInc),wR(posInc),centerGravity,(lenA+lenB),lenA,alphaF,tireCos,camberRear,driveWheels,springF,springR,0,0,aeroCD,aeroCDF,aeroFront,area,rho,1,curveRad);
                   end
                   
                   if(abs(axOld-ax(posInc))<.1)
                       break
                   end
                   
                end

                ay(posInc) = u(posInc).^2./curveRad;
            %% if we aren't in a turn
            else %if we arent in a turn
                    u(posInc) = sqrt(u(posInc-1).^2+2.*ax(posInc-1).*d(posInc));
                    v(posInc) = 0;
                   ay(posInc) = 0;
                   if(driveWheels == 1)  
                        [ax(posInc),wF(posInc),wR(posInc),kappa] = longAccel(u(posInc),v(posInc),gearRatio,finalDrive,torque,wheelRadRear,rpm,wF(posInc),wR(posInc),centerGravity,(lenA+lenB),lenA,0,tireCos,camberRear,driveWheels,springF,springR,0,0,aeroCD,aeroCDF,aeroFront,area,rho,1,NaN);
                   else
                        [ax(posInc),wF(posInc),wR(posInc),kappa] = longAccel(u(posInc),v(posInc),gearRatio,finalDrive,torque,wheelRadRear,rpm,wF(posInc),wR(posInc),centerGravity,(lenA+lenB),lenA,0,tireCos,camberRear,driveWheels,springF,springR,0,0,aeroCD,aeroCDF,aeroFront,area,rho,1,NaN);
                   end
            end
            
            posInc = posInc + 1;
    
        end

end

uBack(apexInd) = u(apexInd);
vBack(apexInd) = v(apexInd);
axBack(apexInd) = ax(apexInd);
ayBack(apexInd) = ay(apexInd);
wFBack(apexInd) = wfs; 
wRBack(apexInd) = wrs; 

%% braking section
for ii = 2:size(apexStart,1) %do this for every apex
   pressIn = pIn;
   posInc = apexStart(ii)-1;
   if (ii ~= 1) %make sure we arent at the start of the track
       endCond = apexEnd(ii-1);
   else
       endCond = 1;
   end
        while(posInc > endCond)
            wFBack(posInc) = wFBack(posInc+1);
            wRBack(posInc) = wRBack(posInc+1);
                        %% if we are in a turn
            if (~isnan(rads(posInc))) 
                %find velocity at new point based on last point data
                uBack(posInc) = sqrt(uBack(posInc+1).^2-2.*axBack(posInc+1).*abs(pos(posInc,1)-pos(posInc+1,1)));
                vBack(posInc) = sqrt(vBack(posInc+1).^2+2.*ayBack(posInc+1).*abs(pos(posInc,2)-pos(posInc+1,2)));
                
                %curve Radius
                curveRad = rads(posInc);
                %find force balence
                SA = 0;

                axBack(posInc) = 1000; %initialize value for use in iterative method
                kappa = 0; %initialize value for use in iterative method
                %% find tire balence
                while(1)
                    axOld = axBack(posInc);
                    
                    [~,~,alphaPeakF,~] = magicTire(wFBack(posInc)/2000,0,0,camberFront,tireCos,1);
                    [~,~,alphaPeakR,~] = magicTire(wRBack(posInc)/2000,0,0,camberRear,tireCos,1);
                   % find the slip angles of the tires so that traction can
                   % be checked
                    fyFront(posInc) = mass.*(uBack(posInc).^2./curveRad)./(1+lenA./lenB);
                    fYF_diff = 2.*magicTire(wFBack(posInc)/2000,alphaSweep,0,camberFront,tireCos,1)-fyFront(posInc);
                    alphaF(posInc) = alphaSweep( find(abs(fYF_diff) == min(abs(fYF_diff))) );
                    if (alphaF(posInc)>alphaPeakF)%if the front tires cant generate the tracion to make the turn
                        alphaF(posInc) = alphaPeakF;
                        fyFront(posInc) = 2.*magicTire(wFBack(posInc)/2000,alphaF(posInc),0,camberFront,tireCos,1);
                        uBack(posInc)  = sqrt(fyFront(posInc).*(1+lenA./lenB).*curveRad./mass);
                    end
                        
                    fyRear(posInc) = lenA./lenB.*fyFront(posInc);
                    fYR_diff = 2.*magicTire(wRBack(posInc)/2000,alphaSweep,kappa,camberRear,tireCos,1)-fyRear(posInc);
                    alphaR(posInc) = alphaSweep( find(abs(fYR_diff) == min(abs(fYR_diff))) );
                   
                   [axBack(posInc),wFBack(posInc),wRBack(posInc),~] = brakingAccel(uBack(posInc),vBack(posInc), pressIn, pressMax, wheelRadFront,wheelRadRear, bgf, bgr, pCritical, proportion, wfs, wrs, centerGravity, (lenA+lenB), lenA, SA, tireCos, camberFront,camberRear, springF,springR, 0, 0, aeroCD, aeroCDF, aeroFront, area, rho, 1, curveRad);
                   
                   if(abs(axOld-axBack(posInc))<.1)
                       break
                   end
                   
                end

                ayBack(posInc-1) = uBack(posInc).^2./curveRad;
            %% if we aren't in a turn
            else %if we arent in a turn
                    uBack(posInc) = sqrt(uBack(posInc+1).^2-2.*axBack(posInc+1).*d(posInc));
                    vBack(posInc) = 0;
                    ayBack(posInc) = 0;
                   
                   [axBack(posInc),wFBack(posInc),wRBack(posInc),~] = brakingAccel(uBack(posInc),vBack(posInc), pressIn, pressMax, wheelRadFront,wheelRadRear, bgf, bgr, pCritical, proportion, wfs, wrs, centerGravity, (lenA+lenB), lenA, 0, tireCos, camberFront,camberRear, springF,springR, 0, 0, aeroCD, aeroCDF, aeroFront, area, rho, 1, NaN);
            end
            
            posInc = posInc - 1;
            if (ii > 1 && posInc == endCond)
                if(uBack(posInc+1)<uBack(apexEnd(ii-1)))
                      if (abs(pressIn-pressMax)>200)
                          pressIn = (pressIn+pressMax)./2;
                          posInc = apexStart(ii)-1;
                          disp('roll back')
                      else
                        endCond = apexEnd(ii-2);
                        disp('extended')
                      end                     
                end
            end
    
        end


end
  
%% merge u profiles
for ii = 1:size(apexEnd,1) %do this for every apex
   posStart = 1 + apexEnd(ii);
   posInc = 1 + apexEnd(ii);
   if (ii ~= size(apexStart,1)) %make sure we arent at the end of the track
       endCond = apexStart(ii+1);
       while(u(posInc)<uBack(posInc))
            posInc = posInc + 1;
        end
       uFinal(posStart:posInc) = u(posStart:posInc);
       uFinal(posInc+1:endCond) = uBack(posInc+1:endCond);
   else
       endCond = size(pos,1);
       uFinal(posInc:endCond) = u(posInc:endCond);
   end

end

%% find time

for ii = 1:size(uFinal,2)-1
   dt(ii) = (2.*d(ii))./(uFinal(ii)+uFinal(ii+1));
   time = time + dt(ii);
   timeAccum(ii) = time;
end

disp(strcat("Time around track = ",num2str(time)," seconds"))

velo = sqrt(u.^2+v.^2);
veloBack = sqrt(uBack.^2+vBack.^2);
% 
% figure(2)
% subplot(1,2,1)
% scatter(1:size(alphaF,2),alphaF)
% subplot(1,2,2)
% scatter(1:size(alphaR,2),alphaR)
% 
% figure(1)
% scatter(1:size(u,2),u)
% hold on
% scatter(1:size(d,2),d)
% 
% figure(3)
% scatter(1:size(fyFront,2),fyFront+fyRear)
% hold on
% scatter(apexInd,ones(length(apexInd),1).*-1000)
% 
% figure(4)
% scatter(1:size(u,2),velo)
% hold on
% scatter(1:size(u,2),u)
% scatter(1:size(u,2),v)
% legend('velo','u','v')
% 
% figure(5)
% subplot(1,2,1)
% scatter(1:size(wF,2),wF)
% title('WF')
% subplot(1,2,2)
% scatter(1:size(wR,2),wR)
% title('WR')
% 
% 
% figure(6)
% subplot(1,2,1)
% scatter(1:size(ay,2),ay)
% title('ay')
% subplot(1,2,2)
% scatter(1:size(ax,2),ax)
% title('a')

figure(1)
scatter(1:size(uBack,2),uBack)
hold on
scatter(1:size(d,2),d)

figure(2)
scatter(1:size(uBack,2),uBack)
hold on
scatter(1:size(u,2),u)

figure(3)
scatter(1:size(veloBack,2),veloBack)
hold on
scatter(1:size(velo,2),velo)

figure(4)
scatter(1:size(uFinal,2),uFinal)

[~,~,uReal] = importRawTrackwSpeed('Log-20191020-134416 winston 8.19.20 - 0.45.949 - Interpolated.csv',11050, 12001);

close all
uFinal = uFinal(1:940);
uReal = uReal(1:940);
trackPer = 1:size(uFinal,2)-1;
trackPer = trackPer./(size(uFinal,2)-1).*100;

figure(1)
scatter(trackPer,uFinal(1:length(uFinal)-1).*2.237)
hold on
scatter(trackPer,uReal(1:length(uReal)-1))
%scatter(apexInd,uFinal(apexInd).*2.237)
legend('sim','real')
ylabel('Speed (MPH)')
xlabel('Track Position (%)')

fitLen = 1:size(pos,1);
fitLen = fitLen';
uFinal=uFinal';