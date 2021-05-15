function [Q_epi_avg,ETsoi_act_avg,Q_surf_avg,STATES,FLUXES,Cont_area,Delta_V]=...
    soil_epikarst_routine(fc,f_red,Vr,Ve,a,Vsoi,Vepi,Kepi,n,P_eff,T_pot,...
    Es_pot,t_wet,Vsoi1_ini,Vsoi2_ini,Vsoi3_ini,Vepi_ini,Conc_flow)

% This function simulates the soil and epikarst routine of the V2Karst
% model and evaluate the daily or sub-daily soil and epikarst water balance 
% for the different model vertical compartments.
%
% USAGE:
% [Q_epi_avg,ETsoi_act_avg,Q_surf_avg,STATES,FLUXES,Cont_area]=...
%     soil_epikarst_routine(fc,f_red,Vr,Ve,a,Vsoi,Vepi,Kepi,n,P_eff,T_pot,...
%     Es_pot,t_wet,Vsoi1_ini,Vsoi2_ini,Vsoi3_ini,Vepi_ini,Conc_flow)
%
% INPUTS
%
% PARAMETERS:
%         fc = vegetation cover fraction [-]                      
%              - if constant                                       - scalar
%              - if variable in time                          - vector(H,1)
%      f_red = reduction factor for transpiration in soil 
%              layer 3                                             - scalar
%          Vr = Maximum storage capacity of the rooting 
%               zone [mm]                                          - scalar
%          Ve = Maximum storage capacity of first soil layer    
%              [mm]                                                - scalar
%           a = Spatial variability parameter [-]                  - scalar
%        Vsoi = Mean soil storage capacity [mm]                    - scalar
%        Vepi = Mean epikarst storage capacity [mm]                - scalar
%        Kepi = Mean epikarst outflow coefficient [d]              - scalar
%           n = number of model vertical compartments [-]          - scalar
%
% INPUT DATA:
%       P_eff = effective precipitation that infiltrates in - vector(H,1)
%               the soil (after accounting for interception
%               and change in snow pack) [mm T-1]                 
%       T_pot = potential transpiration [mm T-1]              - vector(H,1)
%      Es_pot = potential soil evaporation [mm T-1]           - vector(H,1)
%       t_wet = fraction of the time step with wet canopy [-] - vector(H,1)
%
% INITIAL STATES:
%  Vsoi1_ini = initial moisture in soil layer 1 
%              [in % saturation of the deepest soil compartment]   - scalar
%  Vsoi2_ini = initial moisture in soil layer 2                  
%              [in % saturation of the deepest soil compartment]   - scalar 
%  Vsoi3_ini = initial moisture in soil layer 2                     
%              [in % saturation of the deepest soil compartment]   - scalar
%   Vepi_ini = initial moisture in epikarst                         
%              [in % saturation of the deepest soil compartment]   - scalar
% OTHER INPUTS:
%  Conc_flow = flag for activation of the concentration flow       - scalar       
%             component of the model (karst processes)                              
%                 - Conc_flow = 1: compartment saturation excess
%                   is channelled to the next unsaturated 
%                   compartment(concentration flow)
%                 - Conc_run = 0: compartment saturation excess
%                   is lost as runoff (no concentation flow)           
% OUTPUTS
%  Q_epi_avg = average recharge over all compartments [mm T-1]- vector(H,1)              
% ETsoi_act_avg = average sum of actual transpiration         - vector(H,1)
%                 and soil evaporation over all compartments 
%                 [mm T-1]
% Q_surf_avg = average surface runoff over all                - vector(H,1)
%              compartments [mm T-1]
%     STATES = average state variables over all               - vector(H,6)
%              compartments [% saturation]
%              STATES(:,1) = total soil water storage
%              STATES(:,2) = epikarst water storage
%              STATES(:,3) = soil water storage in layers 1+2
%                            (rooting zone)
%              STATES(:,4) = soil water storage in layer 1
%              STATES(:,5) = soil water storage in layer 2
%              STATES(:,6) = soil water storage in layer 3
%     FLUXES = average fluxes over all compartments [mm T-1]  - vector(H,6)              
%              FLUXES(:,1) = actual soil evaporation
%              FLUXES(:,2) = actual transpiration in soil 
%                            layer 1
%              FLUXES(:,3) = actual transpiration in soil 
%                            layer 2
%              FLUXES(:,4) = actual transpiration in soil 
%                            layer 3
%              FLUXES(:,5) = lateral flow (flow from saturated 
%                            to unsaturated compartments)
%              FLUXES(:,6) = soil saturation excess
%  Cont_area = Number of compartments in which                - vector(H,1)
%              the soil generates a saturation excess 
%              flow to the epikarst (contributing areas)[-]
%
% NOTES:
% - Initial condition: all compartments initially store the same depth of 
%   water, except for shallow compartments for which the initial condition 
%   is equal to their storage capacity when the initial condition exceed
%   the storage capacity.
% - Inputs should verify: Ve < Vr < Vsoi_max(n) where Vsoi_max(n) is the 
%   total storage capacity of the deepest soil compartment.
% 
% REFERENCES:
% Hartmann, A., T. Gleeson, R. Rosolem, F. Pianosi, Y. Wada, and T. Wagener 
% (2015), A large-scale simulation model to assess karstic groundwater 
% recharge over Europe and the Mediterranean, Geosci. Model Dev., 8(6), 
% 1729–1746, doi:10.5194/gmd-8-1729-2015.
%
% Penman, H. L. (1950), The dependance of transpiration on weather and soil 
% conditions, J. Soil Sci., 1(1), 74–89, 
% doi:10.1111/j.1365-2389.1950.tb00720.x.
%
% Rimmer, A. and Hartmann, A. (2012), Simplified conceptual structures and 
% analytical solutions for groundwater discharge using reservoir equations, 
% in Water resources management and modeling, edited by P. Nayak , InTech, 
% Kakinada, India., 217–238 2012.
%
% This function is part of the V2Karst model by F. Sarrazin, A. Hartmann, 
% F. Pianosi, R. Rosolem, T. Wagener (2018, Geosci. Model Dev.)
% V2Karst is provided under the terms of the GNU General Public License 
% version 3.0.
% This function was prepared by Fanny Sarrazin (fanny.sarrazin@ufz.de).

%--------------------------------------------------------------------------
% 1. Prepare parameters
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 1.1 Value of soil, epikarst storage capacity and epikarst outflow 
%     coefficient for all vertical compartments
%--------------------------------------------------------------------------
if n~=1 && a~=0 % Heterogeneous subsurface properties

    % Calculation of maximum values of soil and epikarst properties across
    % all compartments
    % Option 1: original formulation in VarKarst (Hartmann et al., 2015)
    % Vsoi_dummy = Vsoi*(a+1);
    % Vepi_dummy = Vepi*(a+1);
    % Kepi_dummy = Kepi_mean*(a+1);
    % Option 2: formulation in V2Karst (Sarrazin et al., 2017)
    Vsoi_dummy = Vsoi*n/sum(((1:n)/n).^a);
    Vepi_dummy = Vepi*n/sum(((1:n)/n).^a);
    Kepi_dummy = Kepi*n/sum(((1:n)/n).^a);
    
    % Soil and epikarst properties for the different compartments
    Vsoi_max = Vsoi_dummy*((1:n)/n).^a; % soil storage capacity [mm]
    Vepi_max = Vepi_dummy*((1:n)/n).^a; % epikarst storage capacity [mm]
    Ke = max(1,Kepi_dummy*((n-(1:n)+1)/n).^a); % epikarst outflow  
                                       % coefficient [d] (higher than 1 d)
    
elseif n==1 || a==0 % Homogeneous subsurface properties
    Ke = Kepi;
    Vsoi_max = Vsoi;
    Vepi_max = Vepi;
    n = 1; % force n equal to 1 to reduce computational time(only one model
           % compartment is needed)
end

%--------------------------------------------------------------------------
% 1.2 Depth of the 3 soil layers 
%--------------------------------------------------------------------------
% Check that Ve < Vr < Vsoi_max(n)
if Ve > Vr || Vr > Vsoi_max(n)
    error('''Ve'' must be lower than ''Vr'' which is turn must be lower than ''Vsoi_max(n)''')
end

% Effective depth of the first soil layer
Vsoi1_max = min(Ve,Vsoi_max,'includenan'); % vector(1,n)
% Effective rooting depth
Vr_eff = min(Vr,Vsoi_max,'includenan'); % vector(1,n)
% Effective depth of soil layer 2
Vsoi2_max=Vr_eff-Vsoi1_max;
% Effective depth of soil below root zone
Vsoi3_max=Vsoi_max-Vsoi2_max-Vsoi1_max;

% Checks for code debugging
% if any(round(Vsoi1_max*10000)<0)  || any(round(Vsoi2_max*10000)<0) || any(round(Vsoi3_max*10000)<0)
%     error('Soil depth is negative')
% end
% if any(round((Vsoi_max-Vsoi1_max-Vsoi2_max-Vsoi3_max)*10000))~=0
%     error('Sum of soil layer depth is not equal to total soil depth')
% end

% Correct possible numerical errors
Vsoi1_max=max(Vsoi1_max,0);
Vsoi2_max=max(Vsoi2_max,0);
Vsoi3_max=max(Vsoi3_max,0);


%--------------------------------------------------------------------------
% 1.3 Prepare vegetated cover fraction
%--------------------------------------------------------------------------
H = size(P_eff,1); % length of simulation horizon
if isscalar(fc);fc = fc*ones(H,1);end

%--------------------------------------------------------------------------
% 2. Initialise variables
%--------------------------------------------------------------------------

Vsoi1       = nan(H,n); % soil water storage in soil layer 1 [mm]
Vsoi2       = nan(H,n); % soil water storage in soil layer 2 [mm]
Vsoi3       = nan(H,n); % soil water storage in soil layer 3 [mm]
Vsoi_perc   = nan(H,1); % soil water storage [% saturation]
Vsoi12_perc = nan(H,1); % soil water storage in root zone [% saturation]
Vsoi1_perc  = nan(H,1); % soil water storage in soil layer 1 [% saturation]
Vsoi2_perc  = nan(H,1); % soil water storage in soil layer 2 [% saturation]
Vsoi3_perc  = nan(H,1); % soil water storage in soil layer 3 [% saturation]
Vepi        = nan(H,n); % epikarst water storage [mm]
Vepi_perc   = nan(H,1); % epikarst water storage [% saturation]
T1_act      = nan(H,n); % actual transpiration in soil layer 1 [mm]
T2_act      = nan(H,n); % actual transpiration in soil layer 2 [mm]
T3_act      = nan(H,n); % actual transpiration in soil layer 3 [mm]
Es_act      = nan(H,n); % actual soil evaporation [mm]
R12         = nan(H,n); % flow from soil layer 1 to layer 2 [mm]
R23         = nan(H,n); % flow from soil layer 2 to layer 3 [mm]
Repi        = nan(H,n); % flow from soil layer 3 to epikarst [mm]
Exc_Epi     = nan(H,n); % epikarst saturation excess [mm]
Q_epi       = nan(H,n); % groundwater recharge [mm]
Q_lat       = zeros(H,n); % lateral flow [mm]
Q_surf_avg  = zeros(H,1); % average surface runoff across all 
                          % compartments [mm]
Delta_V     = nan(H,1); % variation in water storage compared to initial 
                        % condition [mm]
                        
%--------------------------------------------------------------------------
% 3. Soil and epikarst routine
%--------------------------------------------------------------------------

for t=1:H
    
    %----------------------------------------------------------------------
    % 3.1 Initialise state variables
    %---------------------------------------------------------------------
    if t==1 % set initial conditions
        Vsoi1_0 = min(Vsoi1_ini/100*Vsoi1_max(n),Vsoi1_max,'includenan');% soil layer 1
        Vsoi2_0 = min(Vsoi2_ini/100*Vsoi2_max(n),Vsoi2_max,'includenan');% Soil layer 2
        Vsoi3_0 = min(Vsoi3_ini/100*Vsoi3_max(n),Vsoi3_max,'includenan');% Soil layer 3
        Vepi_0  = min(Vepi_ini/100*Vepi_max(n),Vepi_max,'includenan');   % Epikarst
    else % moisture at the end of the previous time step
        Vsoi1_0 = Vsoi1(t-1,:);
        Vsoi2_0 = Vsoi2(t-1,:);
        Vsoi3_0 = Vsoi3(t-1,:);
        Vepi_0  = Vepi(t-1,:);
    end
    
    %----------------------------------------------------------------------
    % 3.2 Soil evaporation:
    %---------------------------------------------------------------------
    
    Es_act_dummy = (1-fc(t))*Es_pot(t)*min(Vsoi1_0./Vsoi1_max,1,'includenan');
    % Evaporation cannot be higher than available moisture Vsoi1_0 + P_eff(t)
    Es_act(t,:) = min(Es_act_dummy,Vsoi1_0 + P_eff(t),'includenan');
    
    % Update available moisture in layer 1 with soil evaporation
    % (soil evaporation has priority over transpiration)
    Vsoi1_dummy = Vsoi1_0 + P_eff(t) - Es_act(t,:);
    
    %----------------------------------------------------------------------
    % 3.3 Transpiration:
    %----------------------------------------------------------------------
    % Rate at which transpiration would occur in layers 1+2 (root zone):
    T12_rate = fc(t).*T_pot(t)*(1 - t_wet(t))*...
        min((Vsoi1_0 + Vsoi2_0)./(Vsoi1_max+Vsoi2_max),1,'includenan');
    
    % Rate at which transpiration would occur in layer 3 ( below root zone) 
    T3_rate = zeros(1,n);
    for i=1:n
        if Vsoi3_max(i)~=0 % if soil is deep enough to have a layer below root zone
            T3_rate(i) = fc(t)*T_pot(t)*(1 - t_wet(t))*...
                min(Vsoi3_0(i)./Vsoi3_max(i),1,'includenan')*f_red; 
            % The rate of transpiration if reduced by the factor f_red to
            % account for the fact that below the root zone moisture is
            % less easily accessed by the roots (Penman, 1950)
        end
    end
    
    % Transpiration occurs where the rate is higher
    T12_rate(T3_rate > T12_rate) = 0;
    T3_rate(T12_rate >= T3_rate) = 0;
    
    % Actual transpiration cannot be higher than available moisture in root zone
    T12_act = min(T12_rate,Vsoi1_dummy + Vsoi2_0,'includenan');
    
    % Transpiration is distributed among layers 1 and 2 proportionally to
    % their moisture depth
    for i=1:n
        if Vsoi1_dummy(i)+Vsoi2_0(i)~=0
            % correct transpiration in soil layer 1 so that it is equal to or lower 
            % than available moisture and total T rate (to avoid numerical
            % errors)
            T1_act(t,i) =...
                min(min(T12_act(i)*Vsoi1_dummy(i)/(Vsoi1_dummy(i)+Vsoi2_0(i)),...
                T12_act(i),'includenan'),Vsoi1_dummy(i),'includenan');
        else
            T1_act(t,i)=0;
        end
    end
    
    % Actual T in soil layer 2
    T2_act(t,:)=min(T12_act - T1_act(t,:),Vsoi2_0,'includenan');
    
    %----------------------------------------------------------------------
    % 3.4 Soil water balance:
    %----------------------------------------------------------------------
    % The flow from one soil layer to next is computed as the layer
    % saturation excess
    
    % Update soil moisture in layer 1
    Vsoi1(t,:) = min(Vsoi1_dummy - T1_act(t,:),Vsoi1_max,'includenan');
    R12(t,:) = max(0,Vsoi1_dummy - T1_act(t,:)-Vsoi1_max,'includenan');
    
    % Update soil moisture in layer 2
    Vsoi2(t,:) = min(Vsoi2_0 + R12(t,:) - T2_act(t,:),Vsoi2_max,'includenan');
    R23(t,:) = max(0,Vsoi2_0 + R12(t,:) - T2_act(t,:)-Vsoi2_max,'includenan');
    
    % Update soil moisture in layer 3
    T3_act(t,:) = min(T3_rate , Vsoi3_0 + R23(t,:),'includenan');
    Vsoi3(t,:) = min(Vsoi3_0 + R23(t,:) - T3_act(t,:),Vsoi3_max,'includenan');
    Repi(t,:) = max(Vsoi3_0 + R23(t,:) - T3_act(t,:) - Vsoi3_max,0,'includenan');
    
    %----------------------------------------------------------------------
    % 3.5 Epikarst water balance:
    % same as in VarKarst (Hartmann et al., 2015)
    %----------------------------------------------------------------------
    Q_epi_dummy = Vepi_0./Ke; % equation for linear reservoir (Rimmer and 
                                % Hartmann, 2012)
    Q_epi(t,:) = min(Q_epi_dummy,Vepi_0 + Repi(t,:),'includenan');
    Vepi(t,:) = min(Vepi_0 + Repi(t,:) - Q_epi(t,:),Vepi_max,'includenan');
    Exc_Epi(t,:) = max(0,Vepi_0 + Repi(t,:) - Q_epi(t,:) - Vepi_max,'includenan');
    
    %----------------------------------------------------------------------
    % 3.6 Generation of lateral flow and runoff:
    % V2Karst implements a corrected version of VarKarst (Hartmann et al., 2015)
    % to avoid numerical errors.
    %----------------------------------------------------------------------
    
    if sum(Exc_Epi(t,:))~=0 
        
        % Index of last compartment with saturation excess
        idx_sat_last = find(Exc_Epi(t,:)>0, 1, 'last' );
        
        if  (idx_sat_last<n && Conc_flow==1) 
            % Activation of concentration flow component
 
            % 1.Distribution of saturation excess among unsaturated
            % compartments: the next unsaturated compartment is filled until
            % saturation is reached, then the subsequent compartment is filled
            % with the remaining excess moisture and so on. If there is still
            % moisture left after all compartments are filled, this moisture
            % excess becomes surface runoff.
            
            % Compartments moisture deficit
            Deficit = Vepi_max+Vsoi_max-Vsoi1(t,:)-Vsoi2(t,:)-Vsoi3(t,:)-Vepi(t,:);
            
            % Total saturation excess
            Excess=sum(Exc_Epi(t,:));
            % Loop over unsaturated compartments
            for i=idx_sat_last+1:n
                % Lateral flow into compartment i
                Q_lat(t,i)=min(Excess,Deficit(i),'includenan');
                Excess = Excess-Q_lat(t,i);
            end
            % The remaining saturation excess becomes surface runoff
            Q_surf_avg(t)=Excess/n; % divided by n to obtain average runoff 
            % over model compartments 
            
            % 2.Update soil and epikarst moisture and correct states so that
            % they do not exceed maximum storage because of numerical errors
            
            % Concentration flow filling soil layer 1
            Q_lat1 = min(Q_lat(t,:),Vsoi1_max-Vsoi1(t,:),'includenan');
            Vsoi1(t,:) = min(Vsoi1(t,:)+Q_lat1,Vsoi1_max,'includenan');
            % Concentration flow filling soil layer 2
            Q_lat2 = min(Q_lat(t,:)-Q_lat1,Vsoi2_max-Vsoi2(t,:),'includenan');
            Vsoi2(t,:) = min(Vsoi2(t,:)+Q_lat2,Vsoi2_max,'includenan');
            % Concentration flow filling soil layer 3
            Q_lat3 = min(Q_lat(t,:)-Q_lat1-Q_lat2,Vsoi3_max-Vsoi3(t,:),'includenan');
            Vsoi3(t,:) = min(Vsoi3(t,:)+Q_lat3,Vsoi3_max,'includenan');
            % Concentration flow filling epikarst
            Q_latEpi = min(Q_lat(t,:)-Q_lat1-Q_lat2-Q_lat3,Vepi_max-Vepi(t,:),'includenan');
            Vepi(t,:) = min(Vepi(t,:)+Q_latEpi,Vepi_max,'includenan');
            
             % Checks for code debugging
%             if Excess<0;error('Excess<0');end
%             
%             if any(round((Q_lat(t,:)-Q_lat1-Q_lat2-Q_lat3-Q_latEpi)*10000)~=0)
%                 error('error in the distribution of concentration flow between soil and epikarst layers')
%             end
            
        elseif  (idx_sat_last==n || Conc_flow==0)
            % all compartments are saturated or concentration flow is not
            % considered: all saturation excess becomes surface runoff
            Q_surf_avg(t)=sum(Exc_Epi(t,:))/n;% divided by n to obtain 
            % average runoff over model compartments
        end
    end
    %----------------------------------------------------------------------
    % 3.7 Assess average moisture  [% saturation]
    %---------------------------------------------------------------------
    Vsoi_perc(t,:) = mean((Vsoi1(t,:)+Vsoi2(t,:)+Vsoi3(t,:))./Vsoi_max,2);
    Vsoi12_perc(t,:) = mean((Vsoi1(t,:)+Vsoi2(t,:))./(Vsoi1_max+Vsoi2_max),2);
    Vsoi1_perc(t,:) = mean(Vsoi1(t,:)./Vsoi1_max,2);
    Vsoi2_perc_dummy = Vsoi2(t,:)./Vsoi2_max;
    Vsoi2_perc_dummy(Vsoi2_max==0)=0;
    Vsoi2_perc(t,:) = mean(Vsoi2_perc_dummy,2);
    Vsoi3_perc_dummy=Vsoi3(t,:)./Vsoi3_max;
    Vsoi3_perc_dummy(Vsoi3_max==0)=0;
    Vsoi3_perc(t,:) = mean(Vsoi3_perc_dummy,2);
    Vepi_perc(t,:) = mean(Vepi(t,:)./Vepi_max,2);
    
    % Checks for code debugging
%     if round((sum(Q_lat(t,:))+Q_surf_avg(t)*n-sum(Exc_Epi(t,:)))*10000)~=0
%         error('error in the runoff routine')
%     end
%     if Conc_flow==0 && any(Q_lat(t,:)~=0);error('Concentation flow is occurring');end

    %----------------------------------------------------------------------
    % 3.8 Assess variation in water storage compared to initial condition
    %---------------------------------------------------------------------
    
    Delta_V(t) = mean(Vsoi1(t,:) + Vsoi2(t,:) + Vsoi3(t,:) + Vepi(t,:) -...
        min(Vsoi1_ini/100*Vsoi1_max(n),Vsoi1_max,'includenan')-...
        min(Vsoi2_ini/100*Vsoi2_max(n),Vsoi2_max,'includenan')- ...
        min(Vsoi3_ini/100*Vsoi3_max(n),Vsoi3_max,'includenan')-...
        min(Vepi_ini/100*Vepi_max(n),Vepi_max,'includenan'));
    
end

%--------------------------------------------------------------------------
% 4. Compute model outputs (average over model compartments)
%--------------------------------------------------------------------------

T1_act_avg = mean(T1_act,2); % Transpiration in layer 1
T2_act_avg = mean(T2_act,2); % Transpiration in layer 2
T3_act_avg = mean(T3_act,2); % Transpiration in layer 3
Es_act_avg = mean(Es_act,2); % Soil evaporation
ETsoi_act_avg = T1_act_avg+T2_act_avg+T3_act_avg+Es_act_avg; % sum of 
% transpiration and soil evaporation

Q_epi_avg  = mean(Q_epi,2); % Recharge

Cont_area=sum(Repi>0,2); % Number of compartments in which the soil 
                         % generates a saturation excess flow to the epikarst

STATES=[Vsoi_perc,Vepi_perc,Vsoi12_perc,Vsoi1_perc,Vsoi2_perc,Vsoi3_perc];
FLUXES=[Es_act_avg,T1_act_avg,T2_act_avg,T3_act_avg,mean(Q_lat,2),mean(Repi,2)];

%--------------------------------------------------------------------------
% 5 Variable checks
%--------------------------------------------------------------------------
if any(isnan((Q_epi_avg)));error('''Q_epi_avg'' contains NaNs');end
if any(isnan((ETsoi_act_avg)));error('''ETsoi_act_avg'' contains NaNs');end
if any(isnan((Q_surf_avg)));error('''Q_surf_avg'' contains NaNs');end
if any(any(isnan(FLUXES)));error('''FLUXES'' contains NaNs');end
if any(any(isnan(STATES)));error('''STATES'' contains NaNs');end
if any(any(isnan(Cont_area)));error('''Cont_area'' contains NaNs');end

if any(any((Q_epi<0)));error('''Q_epi'' contains negative values');end
if any(any((Vsoi1_perc<0)));error('''Vsoi1_perc'' contains negative values');end
if any(any((Vsoi2_perc<0)));error('''Vsoi2_perc'' contains negative values');end
if any(any((Vsoi3_perc<0)));error('''Vsoi3_perc'' contains negative values');end
if any(any((Vepi<0)));error('''Vepi'' contains negative values');end
if any(any((T1_act<0)));error('''T1_act'' contains negative values');end
if any(any((T2_act<0)));error('''T2_act'' contains negative values');end
if any(any((T3_act<0)));error('''T3_act'' contains negative values');end
if any(any((Es_act<0)));error('''Es_act'' contains negative values');end
if any(any((R12<0)));error('''R12'' contains negative values');end
if any(any((R23<0)));error('''R23'' contains negative values');end
if any(any((Repi<0)));error('''Repi'' contains negative values');end
if any(any((Exc_Epi<0)));error('''Exc_Epi'' contains negative values');end
if any(any((Q_lat<0)));error('''Q_lat'' contains negative values');end
if any(any((Q_surf_avg<0)));error('''Q_lat'' contains negative values');end
if any(any(Cont_area<0));error('''Cont_area'' contains negative values');end
