function [ forceOut, stiffnessCo, alphaPeak, kappaPeak ] = magicTire( fz, SA, kappa, camber, co, output )
%This function creates a model for the magic tire formula by Pacejka
%INPUTS:
%   fz - the load at that set of tires in kN
%   slip - angle
%   camber - camber angle in degrees
%   co - the load influenced and camber influece tire coefficents
%   output - 1 = laterial, 3 = corrective moment, and 2 = longitudinal,
%   comb - combined loading or not (boolean)
%   kappa - slip ratio 0-1
%OUTPUTS:
%   forceOut - the lateral, longitudinal, or moment  
%    in (Newtons or Newton-meters)
%


%13 Coefficient Model with combined slip for Fx



%% ------Lateral Force -----%%
d_a = co.a1(1).*fz.^2+co.a2(1).*fz;

e_a = co.a6(1).*fz.^2+co.a7(1).*fz+co.a8(1);

sh_a = co.a9(1).*camber;

sv_a = ( co.a10(1).*fz.^2+co.a11(1).*fz ).*camber;

c_a = 1.3;

bcd_a = co.a3(1).*sin( co.a4(1).*atan(co.a5(1).*fz) );

b_a = bcd_a./(c_a.*d_a);

phi_a = (1-e_a).*(SA+sh_a)+(e_a/b_a).*atan(b_a.*(SA+sh_a));

if kappa ~= 0
    byk = co.byk(2);
    cyk = co.rcy1(2).*fz.^6+ co.rcy2(2).*fz.^5+ co.rcy3(2).*fz.^4+ co.rcy4(2).*fz.^3+ co.rcy5(2).*fz.^2+ co.rcy6(2).*fz+ co.rcy7(2);
    eyk = co.eyk(2);
    shyk = co.shyk(2);
    gyk = @(kappa) cosd(cyk.*atand(byk.*kappa-eyk.*(byk.*kappa-atand(byk.*kappa)))) ./ cosd(cyk.*atand(byk.*shyk-eyk.*(byk.*shyk-atand(byk.*shyk))));
    gyk = gyk(kappa);
else
    gyk = 1;
end

fy = (d_a.*sin(c_a.*atan(b_a.*phi_a))+sv_a).*gyk;

alpha_array = 0:.001:20;
phi_array = (1-e_a).*(alpha_array+sh_a)+(e_a/b_a).*atan(b_a.*(alpha_array+sh_a));
fy_array = d_a.*sin(c_a.*atan(b_a.*phi_array))+sv_a;

alphaPeak = alpha_array(find(max(fy_array)==fy_array));
phi_array = [];
%% -------corrective Moment
d_m = co.a1(3).*fz.^2+co.a2(3).*fz;

e_m = co.a6(3).*fz.^2+co.a7(3).*fz+co.a8(3);

sh_m = co.a9(3).*camber;

sv_m = ( co.a10(3).*fz.^2+co.a11(3).*fz ).*camber;

c_m = 2.4;

bcd_m = (co.a3(3).*fz.^2+co.a4(3).*fz)./(exp(co.a5(3).*fz));

b_m = bcd_m./(c_m.*d_m);

phi_m = (1-e_m).*(SA+sh_m)+(e_m/b_m).*atan(b_m.*(SA+sh_m));

mz = d_m.*sin(c_m.*atan(b_m.*phi_m))+sv_m;



%% ----- longitudinal force
d_k = co.a1(2).*fz.^2+co.a2(2).*fz;

e_k = co.a6(2).*fz.^2+co.a7(2).*fz+co.a8(2);

sh_k = co.a9(2).*camber;

sv_k = ( co.a10(2).*fz.^2+co.a11(2).*fz ).*camber;

kappa = kappa.*100;    
    
c_k = 1.65;

bcd_k = (co.a3(2).*fz.^2+co.a4(2).*fz).*(exp(-co.a5(2).*fz));

b_k = bcd_k./(c_k.*d_k);

phi_k = (1-e_k).*kappa+(e_k/b_k).*atan(b_k.*kappa);

if SA ~= 0
    bxa = co.bxa(2);
    cxa = co.rcx1(2).*fz.^6+ co.rcx2(2).*fz.^5+ co.rcx3(2).*fz.^4+ co.rcx4(2).*fz.^3+ co.rcx5(2).*fz.^2+ co.rcx6(2).*fz+ co.rcx7(2);
    exa = co.exa(2);
    shxa = co.shxa(2);
    gxa = @(alpha) cosd(cxa.*atand(bxa.*alpha-exa.*(bxa.*alpha-atand(bxa.*alpha)))) ./ cosd(cxa.*atand(bxa.*shxa-exa.*(bxa.*shxa-atand(bxa.*shxa))));
    gxa = gxa(SA);
else
    gxa = 1;
end

fx = -d_k.*sin(c_k.*atan(b_k.*phi_k)).*gxa;

kappa_array = 0:.1:50;

phi_array = (1-e_k).*kappa_array+(e_k/b_k).*atan(b_k.*kappa_array);
fx_array = d_k.*sin(c_k.*atan(b_k.*phi_array)).*gxa;

kappaPeak = kappa_array(find(max(abs(fx_array))==abs(fx_array)))./100;



%% -----Force ouput selection ------
if (output == 1) %Laterial Direction
    forceOut = fy;
    stiffnessCo = bcd_a;
elseif( output == 2)%Longitudinal direction
    forceOut = -fx;
    stiffnessCo = bcd_k;
elseif( output == 3)%Corrective moment
    forceOut = mz;
    stiffnessCo = bcd_m; 
end



end

