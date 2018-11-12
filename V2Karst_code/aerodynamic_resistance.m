function ra = aerodynamic_resistance(zom,zoh,d,WS,zm,zh)

% This function calculates the daily or sub-daily aerodynamic resistance 
% following Allen et al., 1998, Eq.(4).
%
% USAGE:
%  ra = aerodynamic_resistance(zom,zoh,d,WS,zm,zh)
%
% INPUTS:
% zom = roughness length for momentum transfer [m]            - scalar
% zoh = roughness length for heat and vapor transfer [m]      - scalar
%   d = zero plane displacement height [m]                    - scalar              
%  WS = daily wind speed [m/s]                                - vector(H,1)
%  zm = height of wind speed measurements [m]                 - scalar
%  zh = height of humidity measurements [m]                   - scalar
% 
% OUTPUTS:
% ra = aerodynamic resistance [s/m]                           - vector(H,1)
%
% REFERENCES:
% Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998),Crop 
% evapotranspiration:  Guidelines for computing crop requirements, FAO 
% Irrigation and Drainage Paper 56, Food and Agriculture Organization 
% (FAO), Rome, Italy
%
% This function is part of the V2Karst model V1.1 by F. Sarrazin, A. 
% Hartmann, F. Pianosi, R. Rosolem, T. Wagener (2019, Geosci. Model Dev.)
% V2Karst is provided under the terms of the GNU General Public License 
% version 3.0.
% This function was prepared by Fanny Sarrazin, University of Bristol,
% November 2018 (fanny.sarrazin@bristol.ac.uk).

% Define constant parameters
k = 0.41; % von Karman constant [-]

% Calculate aerodynamic resistance
ra = log((zm-d)/zom)*log((zh-d)/zoh)./(k^2*WS);

% Check variable
if any(isnan(ra));error('''ra'' contains NaNs');end
