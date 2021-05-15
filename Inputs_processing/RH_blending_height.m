function [RH_b, E_pot]= RH_blending_height(rs,ra,ra_b,Rn,Ta,RH,G,Pa,Kt)

% This function computes rdaily/sub-daily relative humidity at a blending
%  height following Lhomme et al.(2014).
%
% USAGE:
% RH_b = RH_blending_height(rs,ra,ra_b,Rn,Ta,RH,G,Pa)
%
% INPUTS: 
%   rs = surface resistance [s m-1]                 - scalar or vector(H,1)
%   ra = aerodynamic resistance [s m-1]             - scalar or vector(H,1)
% ra_b = aerodynamic resistance at blending         - scalar or vector(H,1)
%        heigth [s m-1]
%   Rn = net radiation [MJ m-2 T-1]                           - vector(H,1)
%   Ta = mean air temperature [deg C]                         - vector(H,1)
%      or matrix of mean (Ta(:,3)),                          or vector(H,3)    
%      minimum(Ta(:,2))and maximum temperature [deg C]   
%     (see Note below)
%   RH = mean relative humidity                               - vector(H,1)
%    G = ground heat flux [MJ m-2 T-1]                        - vector(H,1)
%   Pa = atmospheric pressure [kPa]                 - scalar or vector(H,1)
%   Kt = time step in seconds [s T-1]                        - scalar 
%        for daily time step set to 3600*24 s d-1
% OUTPUTS:
% RH_b = relative humidity at blending height [%]       - vector(H,1)
%
% NOTES ON INPUT DATA: 
% 1. For time time step, it is recommended to compute the saturation vapour 
%    pressure using daily minimum and maximum temperature due to the 
%    non-linearity of the relationships (Allen et al., 1998).
% 2. Ground heat flux can be neglected at daily time scale (Allen et al.,
%    1998; Shuttleworth, 2012).
% 3. RH_b is adjusted so that it is [0, 100].
% 4. The change in atmospheric density (ro_a) between blending height and 
%    measurement height is neglected (we assume the value at measurement height)
%
% REFERENCES:
% Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998),Crop evapotranspiration: 
% Guidelines for computing crop requirements, FAO Irrigation and Drainage 
% Paper 56, Food and Agriculture Organization (FAO), Rome, Italy
%
% Lhomme, J. P., N. Boudhina, and M. M. Masmoudi (2014), Technical Note: 
% On the Matt–Shuttleworth approach to estimate crop water requirements, 
% Hydrol. Earth Syst. Sci., 18(11), 4341–4348, doi:10.5194/hess-18-4341-2014.

%--------------------------------------------------------------------------
% 1. Prepare variables 
%--------------------------------------------------------------------------
nb_Ta=size(Ta,2); % number of colums in temperature input

%--------------------------------------------------------------------------
% 2. Define constant parameters
%--------------------------------------------------------------------------
cp = 1.013*10^(-3); % specific heat of air at constant pressure for moist  
                    % air [MJ kg-1 C-1]
epsilon = 0.622; % ratio molecular weight of water vapour/dry air [-]
R = 0.287; % specific gas constant [kJ/kg/K]

%--------------------------------------------------------------------------
% 3. Assess variables used in Penman Monteith equation
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 3.1 Latent heat of vaporization
%--------------------------------------------------------------------------
% Formula used in Sarrazin et al., (2017) (Shuttleworth, 1993, Eq.(4.2.1))
lambda = 2.501-0.002361*Ta(:,1);% [MJ/kg] 
% Other option: simplication for a temperature of 20 C (Allen et al., 1998, p.3)
% lambda=2.45; % [MJ/kg] 

%--------------------------------------------------------------------------
% 3.2 Psychrometric constant (Shuttleworth, 1993, Eq.(4.2.28); 
%     Allen et al, 1998, Eq.(8))
%--------------------------------------------------------------------------
gamma = cp*Pa./(lambda*epsilon); % [kPa C-1] 

%--------------------------------------------------------------------------
% 3.3 Saturation vapour pressure (Shuttleworth, 1993, Eq.(4.2.2); 
%     Allen et al., 1998, Eq. (11-12))
%--------------------------------------------------------------------------
if nb_Ta==1
    es = 0.6108*exp((17.27*Ta(:,1))./(Ta(:,1)+237.3));% [kPa] 
elseif nb_Ta==3
    es_min = 0.6108*exp((17.27*Ta(:,2))./(Ta(:,2)+237.3)); % [kPa]
    es_max = 0.6108*exp((17.27*Ta(:,3))./(Ta(:,3)+237.3)); % [kPa]
    es = (es_min+es_max)/2; % [kPa]
end

%--------------------------------------------------------------------------
% 3.4 Actual vapor pressure (Allen et al., 1998, Eq.(17-19))
%--------------------------------------------------------------------------
ea = RH.*es/100; % [kPa]

%--------------------------------------------------------------------------
% 3.5 Gradient of the saturated vapour pressure-temperature function 
% (Shuttleworth, 1993, Eq.(4.2.3); Allen et al., 1998, Eq.(13))
%--------------------------------------------------------------------------
delta = 4098*0.6108*exp((17.27*Ta(:,1))./(Ta(:,1)+237.3))./((Ta(:,1)+237.3).^2); 
% [kPa C-1] 

%--------------------------------------------------------------------------
% 3.6 Atmospheric density (Allen et al., 1998, Eq.(3.6))
%--------------------------------------------------------------------------
TK = (Ta(:,1)+273.15);% temperature [K], uses mean air temperature
TKv = TK./(1-0.378*ea./Pa); % Virtual temperature [K]
ro_a = Pa./(R*TKv); % atmospheric density [kg m-3]

%--------------------------------------------------------------------------
% 4. Compute RH at blending heigth (Lhomme et al., 2014, Eq. (5))
%--------------------------------------------------------------------------
fact = ((delta+gamma).*ra_b+gamma.*rs)./((delta+gamma).*ra+gamma.*rs);
% Vapour pressure deficit at blending height:
D_b = ((es-ea)+delta.*(Rn-G).*ra./(ro_a.*cp*Kt)).*fact-delta.*(Rn-G).*ra_b./(ro_a.*cp*Kt);
% Relative humidity at blending height:
RH_b = (1-D_b./es)*100; 

%--------------------------------------------------------------------------
% 4. Compute potential evapotranspiration 
%--------------------------------------------------------------------------
% Potential evapotranspiration at measurement height:
E_pot = (delta.*(Rn-G)+Kt*ro_a.*cp.*(es-ea)./ra)./(lambda.*(delta+gamma.*(1+rs./ra)));
% Potential evapotranspiration at blending height:
ea_b = RH_b.*es/100;
E_pot_b = (delta.*(Rn-G)+Kt*ro_a.*cp.*(es-ea_b)./ra_b)./(lambda.*(delta+gamma.*(1+rs./ra_b)));

% Check that potential evapotranspiration at measurement height and at
% blending height are equal:
if max(abs(E_pot-E_pot_b))>10^(-5)
    error('Potential evapotranspiration is not the same at measurement height and blending height')
end

% Set to 0 potential evapotranspiration when it is negative:
E_pot = max(E_pot, 0);

%--------------------------------------------------------------------------
% 5. Check variables
%--------------------------------------------------------------------------
% Correct RH_b so that it is between 0 and 100:
if any(RH_b<0)
%     warning('%d values of ''RH_b'' were <0 and set to 0',sum(RH_b<0))
    RH_b(RH_b<0)=0;
end
if any(RH_b>100)
 %   warning('%d values of ''RH_b'' were >100 and set to 100',sum(RH_b>100))
    RH_b(RH_b>100)=100;
end
if any(isnan(RH_b));error('RH_b is NaN');end

% if any(lambda<0); error('lambda is negative');end
% if any(gamma<0); error('gamma is negative');end
% if any(es<0); error('es is negative');end
% if any(ea<0); error('ea is negative');end
% if any(round(ea*10^3)>round(es*10^3)); error([ num2str(sum(ea>es)),' ea must be below es']);end
% if any(delta<0); error('delta is negative');end
% if any(ro_a<0); error('ro_a is negative');end
% if any(isnan(lambda)); error('lambda is NaN');end
% if any(isnan(gamma)); error('gamma is NaN');end
% if any(isnan(es)); error('es is NaN');end
% if any(isnan(ea)); error('ea is NaN');end
% if any(isnan(delta)); error('delta is NaN');end
% if any(isnan(ro_a)); error('ro_a is NaN');end
