function E_pot = Penman_Monteith(rs,ra,Rn,Ta,RH,G,Pa,Kt)

% This function calculates potential evapotranspiration rate 
% using the Penman Monteith equation (Monteith, 1965) for daily or subdaily
% time step.
%
% USAGE: 
% E_pot = Penman_Monteith(rs,ra,Rn,Ta,RH,G,Pa,Kt)
%
% INPUTS: 
% rs = surface resistance [s/m]                     - scalar or vector(H,1)
% ra = aerodynamic resistance [s/m]                 - scalar or vector(H,1)
% Rn = daily net radiation [MJ/m2/step]                       - vector(H,1)
% Ta = mean air temperature                                   - vector(H,1)
%      or matrix of mean (Ta(:,1)),                          or vector(H,3)    
%      minimum(Ta(:,2))and maximum (Ta(:,3)temperature [°C]   
%     (see Note below)
% RH = mean relative humidity                                 - vector(H,1)
%      or matrix of mean (RH(:,1)),                          or vector(H,3) 
%      minimum (RH(:,2))and maximum(RH(:,3)) 
%      relative humidity [%]   
%      (see Note below)
%  G = daily ground heat flux [MJ/m2/step]                    - vector(H,1)
% Pa = atmospheric pressure [kPa]                   - scalar or vector(H,1)
% Kt = conversion factor to estimate aerodynamic                    -scalar 
%      component (Allen et al., 1998 p.26) [s/step]
%      At daily time step 3600*24 s/d
%
% OUTPUTS:
% E_pot = daily potential evapotranspiration [mm]               - vect(H,1)
%
% NOTES ON INPUT DATA: 
% - It is recommended to compute the saturation vapour pressure using daily 
% minimum and maximum temperature and actual vapour pressure using daily 
% minimum and maximum relative humidity due to the non-linearity of the 
% relationships (Allen et al., 1998).
% - Ground heat flux can be neglected at daily time scale (Allen et al.,
% 1998; Shuttleworth, 2012).
%
% REFERENCES:
% Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998),Crop evapotranspiration: 
% Guidelines for computing crop requirements, FAO Irrigation and Drainage 
% Paper 56, Food and Agriculture Organization (FAO), Rome, Italy.
%
% Monteith, J. L.(1965), Evaporation and environment, Symp. Soc. Exp. Biol., 
% 19, 205–234.
%
% Sarrazin, F., A. Hartmann, F. Pianosi, and T. Wagener (2017), V2Karst: 
% A parsimonious large-scale integrated vegetation-recharge model to  
% simulate the impact of climate and land cover change in karst regions. 
% Geosci. Model Dev. In review.
%
% Shuttleworth, W. J.(1993), Evapotranspiration, in Handbook of Hydrology, 
% edited by D. R. Maidment, p. 4.1-4.53, McGraw-Hill inc., New York.
%
% Shuttleworth, J. W. (2012), Terrestrial Hydrometeorology, John Wiley & 
% Sons, Ltd, Chichester, UK.

% This function is part the V2Karst model (Sarrazin et al., 2017). 
% V2Karst is provided under the terms of the GNU General Public License 
% version 3.0.
% This function was prepared by Fanny Sarrazin, University of Bristol,
% August 2018 (fanny.sarrazin@bristol.ac.uk).

%--------------------------------------------------------------------------
% 1. Prepare variables 
%--------------------------------------------------------------------------
nb_Ta=size(Ta,2); % number of colums in temperature input
nb_RH=size(RH,2); % number of colums in relative humidity input

%--------------------------------------------------------------------------
% 2. Define constant parameters
%--------------------------------------------------------------------------
cp = 1.013*10^(-3); % specific heat of air at constant pressure for moist  
                    % air [MJ/kg/°C]
epsilon = 0.622; % ratio molecular weight of water vapour/dry air [-]
R = 0.287; % specific gas constant [kJ/kg/K]

%--------------------------------------------------------------------------
% 3. Assess variables used in Penman Monteith equation
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 3.1 Latent heat of vaporization
%--------------------------------------------------------------------------
% Formula used in Sarrazin et al., (2018) (Shuttleworth, 1993, Eq.(4.2.1))
lambda = 2.501-0.002361*Ta(:,1);% [MJ/kg] 
% Other option: simplication for a temperature of 20°C (Allen et al., 1998, p.3)
% lambda=2.45; % [MJ/kg] 

%--------------------------------------------------------------------------
% 3.2 Psychrometric constant (Shuttleworth, 1993, Eq.(4.2.28); 
%     Allen et al, 1998, Eq.(8))
%--------------------------------------------------------------------------
gamma = cp*Pa./(lambda*epsilon); % [kPa/°C] 

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
if nb_RH==1 || nb_Ta==1
    ea = RH(:,1).*es/100; % [kPa]
elseif nb_RH==3 && nb_Ta==3 % [kPa]
    ea = (RH(:,2).*es_max+RH(:,3).*es_min)/100/2;
end

%--------------------------------------------------------------------------
% 3.5 Gradient of the saturated vapour pressure-temperature function 
% (Shuttleworth, 1993, Eq.(4.2.3); Allen et al., 1998, Eq.(13))
%--------------------------------------------------------------------------
delta = 4098*0.6108*exp((17.27*Ta(:,1))./(Ta(:,1)+237.3))./((Ta(:,1)+237.3).^2); 
% [kPa/°C] 

%--------------------------------------------------------------------------
% 3.6 Atmospheric density (Allen et al., 1998, Eq.(3.6))
%--------------------------------------------------------------------------
TK = (Ta(:,1)+273.15);% temperature [K], uses mean air temperature
TKv = TK./(1-0.378*ea./Pa); % Virtual temperature [K]
ro_a = Pa./(R*TKv); % atmospheric density [kg/m3]

%--------------------------------------------------------------------------
% 4. Assess Potential evapotranspiration
%--------------------------------------------------------------------------
E_pot = max((delta.*(Rn-G)+Kt*ro_a.*cp.*(es-ea)./ra)./(lambda.*(delta+gamma.*(1+rs./ra))),0);

% Check variables
if any(isnan(E_pot));error('''E_pot'' contains NaNs');end
