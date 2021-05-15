function Rn = net_radiation(alpha,emissivity,SW_down,LW_down,Ta)

% This function computes net radiation from downward longwave and shortwave
% radiation, surface albedo and air temperature as a proxy for surface
% temperature.
% 
% USAGE:
% Rn = Net_radiation(alpha,SW_down,LW_down,Ta)
%
% INPUTS:
%      alpha = surface albedo [-]                                
%              - constant                                     - scalar
%              - OR variable in time                          - vector(H,1)
% emissivity = surface emissivity [-]                         
%              - constant                                     - scalar
%              - OR variable in time                          - vector(H,1)
%    SW_down = short wave downward radiation                  - vector(H,1)
%              [MJ m-2 T-1]
%    LW_down = long wave downward radiation                   - vector(H,1)
%              [MJ m-2 T-1]
%          Ta = mean air temperature [deg C]                  - vector(H,1)  
%
% OUTPUTS:
%      Rn = net radiation [MJ m-2 T-1]                        - vector(H,1)
%
% NOTE: The assumption that air temperature is a proxy for surface
%       temperature has only been tested for daily time step. We recommend
%       to use this function for daily time step only and to use a
%       different strategy for subdaily time step.
%
% This function is part of the V2Karst model by F. Sarrazin, A. Hartmann, 
% F. Pianosi, R. Rosolem, T. Wagener (2018, Geosci. Model Dev.)
% V2Karst is provided under the terms of the GNU General Public License 
% version 3.0.
% This function was prepared by Fanny Sarrazin (fanny.sarrazin@ufz.de).

%--------------------------------------------------------------------------
% 1. Check inputs
%--------------------------------------------------------------------------
[H,~]=size(Ta); 
% if  size(SW_down,2)~=1 || size(LW_down,2)~=1  || size(Ta,2)~=1 
%     error('''SW_down'',''LW_down'' and ''Ta'' must be column vectors')
% end
% if size(LW_down,1)~=H || size(SW_down,1)~=H 
%     error('''SW_down'',''LW_down'' and ''Ta'' must have the same number of rows')
% end

if any(alpha)<0 || any(alpha)>1
    error('''alpha'' must be between 0 and 1')
end
if isscalar(alpha);alpha = alpha*ones(H,1);end
if any(emissivity)<0 || any(emissivity)>1
    error('''alpha'' must be between 0 and 1')
end
if isscalar(emissivity);emissivity = emissivity*ones(H,1);end

%--------------------------------------------------------------------------
% 2 Define constants
%--------------------------------------------------------------------------
sigma = 5.67*10^(-8)*3600*24/10^6; % [MJ.m-2.K-4.d-1] Stefan Boltzmann constant 

%--------------------------------------------------------------------------
% 3. Assess net radiation
%--------------------------------------------------------------------------
% Net short wave radiation 
Rns = SW_down.*(1-alpha); % [MJ m-2 T-1] Net shortwave radiation

% Net long wave radiation 
% Net longwave radiation is derived using the Stefan-Boltzmann law assuming
% that the surface temperature is equal to the air temperature. This
% equation was tested at FLUXNET sites and appeared to be reasonable at
% daily time scale (Sarrazin et al.,in prep.)
% Option 1: Equation using mean air temperature
Rnl = LW_down-emissivity*sigma.*(Ta+273.15).^4; % [MJ m-2 T-1] Net longwave radiation

% Option 2: An formulation using minimum and maximum daily air temperature 
% could also be considered but we did not find any significant improvement 
% at FLUXNET sites compared to option 1(where Ta(:,1) is the minimum temperature and 
% Ta(:,2) is the maximum temperature).
% Rnl = LW_down-emissivity*sigma.*((Ta(:,1)+273.15).^4+(Ta(:,2)+273.15).^4)/2; % [MJ m-2 T-1] Net longwave radiation

% Net radiation
Rn = Rns + Rnl; % [MJ m-2 T-1] Net radiation

%--------------------------------------------------------------------------
% 3. Check net radiation
%--------------------------------------------------------------------------
if any(isnan(Rn)); error('Rn contains NaNs');end
if any(Rn<-500);error('Rn<-500');end
if any(Rn>500);error('Rn>500');end

