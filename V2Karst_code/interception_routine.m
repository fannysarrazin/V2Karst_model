function [Tf,t_wet,Ecan_act] = interception_routine(Vcan,LAI,fc,P,Ecan_pot)

% The function estimates daily actual evaporation from canopy interception.
% The formulation is similar to the daily threshold model proposed by 
% (Savenije, 1997; de Groen, 2002; Gerrits, 2010) and do not consider the
% carry-over of interception from one day to the next (all precipitation 
% which is not evaporated reaches the ground as throughfall).
%
% USAGE:
% [Tf,t_wet,Ecan_act] = interception_routine(Vcan,LAI,fc,P,Ecan_pot)
%
% INPUTS
% PARAMETERS:
%       LAI = Canopy leaf area index [m2/m2]
%             - if constant                                   - scalar
%             - if variable in time                           - vector(H,1)
%      Vcan = Canopy storage capacity [mm per unit of LAI]    - scalar
%        fc = vegetation cover fraction [-]                      
%             - if constant                                   - scalar
%             - if variable in time                           - vector(H,1)
%
% INPUT DATA:
%         P = daily precipitation                             - vector(H,1)
%  Ecan_pot = daily potential evaporation from canopy     
%             interception [mm]                               - vector(H,1)
%
% OUTPUTS:
%       Tf = daily throughfall (fraction of precipitation 
%            that reaches the ground) [mm]                    - vector(H,1)                   
%    t_wet = fraction of the day with wet canopy [-]          - vector(H,1)
%            (this will be used to assess transpiration
%             in the funciton soil_epikarst_routine.m)
% Ecan_act = daily actual evaporation from canopy
%             interception [mm]                               - vector(H,1)
%
% REFERENCES:
% Bohn, T. J., and E. R. Vivoni (2016), Process-based characterization of 
% evapotranspiration sources over the North American monsoon region, Water 
% Resour. Res., 52(1), 358�384, doi:10.1002/2015WR017934.
%
% De Groen, M. M.(2002), Modelling interception and transpiration at monthly 
% time steps: introducing daily variability through Markov chains, Ph.D. 
% thesis, Delft University of Technology, Delft, The Netherlands.
%
% Gerrits, M.(2010), The role of interception in the hydrological cycle, Ph.D. 
% thesis, Delft University of Technology, Delft, The Netherlands.
% 
% Kergoat, L. (1998), A model for hydrological equilibrium of leaf area 
% index on a global scale, J. Hydrol.,212�213, 268�286, 
% doi:10.1016/S0022-1694(98)00211-X.
%
% Sarrazin, F., A. Hartmann, F. Pianosi, and T. Wagener (2017), V2Karst: 
% A parsimonious large-scale integrated vegetation-recharge model to  
% simulate the impact of climate and land cover change in karst regions. 
% Geosci. Model Dev. In review.
%
% Savenije, H. H. G. (1997), Determination of evaporation from a catchment 
% water balance at a monthly time scale, Hydrol. Earth Syst. Sci., 1(1), 
% 93�100, doi:10.5194/hess-1-93-1997.

% This function is part the V2Karst model (Sarrazin et al., 2017). 
% V2Karst is provided under the terms of the GNU General Public License 
% version 3.0.
% This function was prepared by Fanny Sarrazin, University of Bristol,
% December 2017 (fanny.sarrazin@bristol.ac.uk).

%--------------------------------------------------------------------------
% 1. Prepare variables
%--------------------------------------------------------------------------
H = size(P,1);
if isscalar(LAI);LAI = LAI*ones(H,1);end
if isscalar(fc);fc = fc*ones(H,1);end

%--------------------------------------------------------------------------
% 2. Initialise variables
%--------------------------------------------------------------------------

Ecan_act = nan(H,1); % Actual evaporation from interception [mm]
Tf = nan(H,1); % Throughfall [mm]
t_wet = nan(H,1); % Fraction of day with wet canopy [-]

%--------------------------------------------------------------------------
% 3. Interception routine
%--------------------------------------------------------------------------

for t=1:H
    % Effective canopy storage capacity over the vegetated fraction
    %(cell average LAI is rescaled to the vegetated fraction following 
    % Bohn and Vivoni, 2016)
    Vcan_max = Vcan*(LAI(t)./fc(t)); % [m2/m2] over vegetated fraction
    
    % Evaporation from interception store is limited by the potential rate
    % and storage capacity
    Ecan_act_dummy = min(Ecan_pot(t),Vcan_max,'includenan');
    
    % Correction (Ecan_act cannot be higher than interception precipitation) 
    % and area-weighting by fc since interception occurs over the vegetated
    % fraction only
    Ecan_act(t) = fc(t).*min(Ecan_act_dummy,P(t),'includenan');
    
    % Throughfall
    Tf(t) = P(t)-Ecan_act(t);
  
    % Assess the fraction of day with wet canopy 
    if Ecan_pot(t)==0 
        t_wet(t) = 0; % the canopy is assumed dry during the entire day
    else
        t_wet(t) = Ecan_act(t)/(Ecan_pot(t).*fc(t)); % similar to Kergoat (1998)
    end
end

%--------------------------------------------------------------------------
% 4. Check variables
%--------------------------------------------------------------------------
if any(isnan(t_wet));error('t_wet contains NaNs');end
if any(isnan(Tf));error('Tf contain NaNs');end
if any(isnan(Ecan_act));error('Ecan_act contains NaNs');end

if any(t_wet>1 | t_wet<0);error('values of t_wet are not between 0 and 1');end
if any(Tf<0);error('Tf contains negative values');end
if any(Ecan_act<0);error('Ecan_act contains negative values');end
if any(Tf>P);error('Tf contains values higher than P');end
