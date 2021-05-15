function [Pre_eff,T_pot_eff, Es_pot_eff,Swe_new] = snow_routine(Swe,Tf,T_pot,Es_pot,swe_thres)

% This function assesses effective precipitation (Pre_eff) from the variation
% of the snow pack and correct T_pot_eff and Es_pot_eff depending on the 
% snowpack. Potential transpiration and potential soil evaporation are set 
% to 0 in presence of snow pack.
%
% USAGE: 
% [Pre_eff,T_pot_eff, Es_pot_eff,Swe_new] = snow_routine(Swe,Tf,T_pot,Es_pot,swe_thres)
%
% INPUTS:
%        Swe = daily snow water equivalent [mm]               - vector(H,1)
%         Tf = daily throughfall (fraction of precipitation 
%             that reaches the ground) [mm d-1]               - vector(H,1) 
%      T_pot = potential transpiration [mm d-1]               - vector(H,1)
%     Es_pot = potential soil evaporation [mm d-1]            - vector(H,1)
%  swe_thres = threshold on swe above which we consider the   - scalar
%              presence of the snow pack and set the potential 
%              transpiration and soil evaporation to 0 [mm]
% 
% OUTPUTS:
%   Pre_eff = effective precipitation that infiltrates in the - vector(H,1)
%             soil [mm T-1]
% T_pot_eff = effective potential transpiration corrected     - vector(H,1)
%             accounting for the presence of the snow pack 
%             [mm T-1]
% Es_pot_eff = effective potential soil evaporation corrected - vector(H,1)
%             accounting for the presence of the snow pack 
%             [mm T-1]
% Swe_new = effective snow pack corrected for consistency     - vector(H,1)
%           with V2Karst simulations
%
% This function is part of the V2Karst model by F. Sarrazin, A. Hartmann, 
% F. Pianosi, R. Rosolem, T. Wagener (2018, Geosci. Model Dev.)
% V2Karst is provided under the terms of the GNU General Public License 
% version 3.0.
% This function was prepared by Fanny Sarrazin (fanny.sarrazin@ufz.de).

%--------------------------------------------------------------------------
% 1. Prepare variables 
%--------------------------------------------------------------------------
H = size(Swe,1);

%--------------------------------------------------------------------------
% 2. Compute the variation of snowpack from one time step to the next
%--------------------------------------------------------------------------
delta_swe_dummy = [0;Swe(2:end)-Swe(1:end-1)]; % initial variation of snow
% pack is set to 0

% Correct delta_swe since the increase in the snow pack cannot be higher
% than Tf (delta_swe<Tf). The residual increase in snow pack is reported 
% to the next day.
% Report the increase in snowpack that exceeds Tf (delta_swe>Tf)to next daya 
delta_swe = zeros(H,1);
for t=2:H
    delta_swe(t)=delta_swe_dummy(t)+max(delta_swe(t-1)-Tf(t-1),0,'includenan');
end
% Correct delta_swe (so that delta_swe<Tf)
delta_swe = min(delta_swe,Tf,'includenan');

%--------------------------------------------------------------------------
% 3. Correct the snow pack so that is it consistent with the simulation
% results
%--------------------------------------------------------------------------
% Correct snow pack
Swe_new = Swe(1) + cumsum(delta_swe);

%--------------------------------------------------------------------------
% 4. Compute the effective precipitation (Precipitation actually reaching
% the ground)
%--------------------------------------------------------------------------
 Pre_eff = Tf - delta_swe;

%--------------------------------------------------------------------------
% 5. Correct potential evapotranspiration
%--------------------------------------------------------------------------
% No transpiration or soil evaporation occur in presence of snow pack
% (potential transpiration and potential soil evaporation are set to 0 in
% presence of snow pack)
T_pot_eff = T_pot;
Es_pot_eff = Es_pot;
T_pot_eff(Swe_new>swe_thres) = 0;
Es_pot_eff(Swe_new>swe_thres) = 0;

%--------------------------------------------------------------------------
% 6. Check variables
%--------------------------------------------------------------------------
% Check that Pre_eff is positive and not NAN
if any(Pre_eff<0);error('Pre_eff<0');end
if any(isnan(Pre_eff));error('Pre_eff is nan');end
% Check that no snow is lost
if round((sum(Tf-Pre_eff)-(Swe_new(end)-Swe_new(1)))*10^4)/10^4~=0;error('Error in snow routine');end
