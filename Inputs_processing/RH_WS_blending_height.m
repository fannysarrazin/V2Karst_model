function [RH_b,WS_b,Rn,E_pot_ref] = RH_WS_blending_height(Pa,Rad,Ta,RH,G,WS,zm,zh,zb,Kt)

% This function computes relative daily/sub-daily humidity and wind speed 
% at a blending height assuming that weather variables were measured above 
% a reference grass with albedo 0.23, height 0.12 m and surface resistance 
% 70 s/m.
% 
% USAGE:
% [RH_b,WS_b,Rn,E_pot_ref] = RH_WS_blending_height(Pa,Rad,Ta,RH,G,WS,zm,zh,zb,Kt)
%
% INPUTS:
%      Pa = atmospheric pressure[kPa]                     - vector(H,1)
%           (see Note 1 below)
%     Rad = radiation input [MJ m-2 T-1]                      - vector(H,1)
%           Net radiation is either provided as input         - vector(H,2)
%           or calculated based on downward radiation.
%           - option 1: vector of size(H,1) that       
%             contains values of net radiation
%           - option 2: vector of size(H,2) that contains 
%              values of downward radiation. With this 
%              option, net radiation will be computed
%              based on albedo and air temperature
%              as a proxy for surface temperature 
%              (see help of net_radiation.m for further 
%               explanation).
%               *Rad(:,1): short wave downward radiation
%               *Rad(:,2): long wave downward radiation
%            (see Note 2 below)
%      Ta = mean air temperature                              - vector(H,1)
%           or matrix of mean (Ta(:,1)),                     or vector(H,3)    
%           minimum(Ta(:,2))and maximum (Ta(:,3)temperature [C]   
%           (see Note 3 below)
%      RH = mean relative humidity [%]                        - vector(H,1)
%       G = ground heat flux [MJ m-2 T-1]                     - vector(H,1)
%      WS = wind speed at 10 m height [m s-1]                 - vector(H,1)
%      zm = height of wind speed measurements [m]             - scalar
%      zh = height of humidity measurements [m]               - scalar
%      zb = blending height at which relative humidity and    - scalar
%           wind speed must be estimated [m]
%      Kt = time step in seconds [s T-1]                      - scalar 
%           for daily time step set to 3600*24 s d-1
% OUTPUTS:
%    RH_b = mean relative humidity at blending height [%]     - vector(H,1)  
%           (see Note 4 below)
%    WS_b = wind speed at blending heigth [m s-1]             - vector(H,1)
%      Rn = net radiation [MJ m-2 T-1]                        - vector(H,1)
%Epot_ref = potential evapotranspiration for reference        - vector(H,1)
%           grass [mm T-1] 
%
% NOTE:
% 1. When Pa is not available, it can be estimated based on elevation as 
%    reported in Allen et al., 1998, Eq.(7) (this calculation has to be 
%    performed prior to calling RH_WS_blending_height).
% 2. Net radiation in assessed from downward radiations assuming that air 
%    temperature is a proxy for surface temperature, which has only been 
%    tested for daily time step. When net radiation has to be calculated
%    from downward radiations, we recommend to use a daily time step (and
%    NOT sub-daily) to simulate V2Karst.
% 3. It is recommended to assess the Penman Monteith potential 
%    evapotranspiration using minimum and maximum daily temperature due to 
%    the non-linearity of the relationships (Allen et al., 1998).
% 4. RH_b is adjusted so that it is [0, 100]. The change in atmospheric 
%    density (ro_a) between blending height and measurement height is
%    neglected (we assume the value at measurement height).
% 
% REFERENCES:
% Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998),Crop evapotranspiration: 
% Guidelines for computing crop requirements, FAO Irrigation and Drainage 
% Paper 56, Food and Agriculture Organization (FAO), Rome, Italy.
%
% Lhomme, J. P., N. Boudhina, and M. M. Masmoudi (2014), Technical Note: 
% On the Matt–Shuttleworth approach to estimate crop water requirements, 
% Hydrol. Earth Syst. Sci., 18(11), 4341–4348, doi:10.5194/hess-18-4341-2014.

%--------------------------------------------------------------------------
% 1. Check inputs
%--------------------------------------------------------------------------
if  size(Pa,2)~=1 || size(G,2)~=1 || size(WS,2)~=1 || size(RH,2)~=1
    error('''Pa'', ''G'' ''RH'' and''WS'' must be column vectors')
end
if size(Rad,2)>2
    error('''Rad'' must have 1 or 2 columns')
end
if ~any(size(Ta,2)==[1 3])
    error('''Ta'' must have 1 or 3 columns')
end

H=size(WS,1); % length of simulation period
if size(Rad,1)~=H || size(G,1)~=H || size(WS,1)~=H || size(Ta,1)~=H || ...
        size(RH,1)~=H || size(Pa,1)~=H
    error('''Rad'',''Pa'',''G'',''WS'',''Ta''and ''RH'' must have the same number of rows')
end

if ~isscalar(zm); error('''zm'' must be scalar and positive');end
if ~isscalar(zh); error('''zh'' must be scalar and positive');end
if ~isscalar(zb); error('''zb'' must be scalar and positive');end

%--------------------------------------------------------------------------
% 2 Define properties of reference grass (Allen et al., 1998)
%--------------------------------------------------------------------------
albedo = 0.23; % [-] surface albedo for reference grass 
h_veg = 0.12; % [m] reference crop height
rs = 70; % [s/m] reference crop surface resistance
emissivity = 1; % [-] surface emissivity 

if h_veg>zm || h_veg>zh || h_veg>zb; error('zm, zh and zb must be higher than h_veg=0.12');end

%--------------------------------------------------------------------------
% 3. Assess net radiation
%--------------------------------------------------------------------------
if size(Rad,2)==2
    % Net radiation is assessed from downward radiation using the average 
    % albedo over the vegetated and non-vegetated fracions. Use mean temperature
    Rn = net_radiation(albedo,emissivity,Rad(:,1),Rad(:,2),Ta(:,1));
else
    % Net radiation is directly provided as input
    Rn = Rad;
end

%--------------------------------------------------------------------------
% 4. Wind speed at blending height and aerodynamic resistances
%--------------------------------------------------------------------------
% Canopy roughness lengths used in Sarrazin et al. (2017)
%(Allen et al., 1998, BOX 4, p.21)
zom = 0.123*h_veg; % canopy roughness length for momentum transfer [m]
zoh = 0.1*zom; % canopy roughness length for heat and vapor transfer [m] 

% Other option to assess canopy roughness lengths
% (Neitsch et al., 2009, p128)
% if h_veg > 2
%     zom_can=(0.058*(h_veg*100)^1.19)/100; 
% else
%     zom_can = 0.123*h_veg;
% end
% zoh_can = 0.1*zom_can;

% Canopy zero plane displacement height(Allen et al., 1998, BOX 4, p.21)
d = 2/3*h_veg; %  [m]  

% Wind speed at blending height
WS_b = WS_blending_height(zom,d,WS,zm,zb); % [s/m] wind speed

% Canopy aerodynamic resistance [s/m]
ra = aerodynamic_resistance(zom,zoh,d,WS,zm,zh); % [s/m] at measurement height
ra_b = aerodynamic_resistance(zom,zoh,d,WS_b,zb,zb); % [s/m] at blending height

%--------------------------------------------------------------------------
% 6. Relative humidity at blending height and potential evapotranspiration
%    for reference grass
%--------------------------------------------------------------------------
[RH_b, E_pot_ref] = RH_blending_height(rs,ra,ra_b,Rn,Ta,RH,G,Pa,Kt); 
 
