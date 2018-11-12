function [Q_epi_d,ET_act_d,Q_surf_d,FLUXES_d,STATES_d,...
          Cont_area_d,time_sim] =...
          V2Karst(param,n,zm,zh,z,P,Rn,Ta,RH,G,WS,time,...
          States_ini,Conc_flow,warmup,Kt)

% This function simulates the V2Karst model V1.1, which is a vegetation-
% recharge model for karst areas. The model solves the water balance at 
% daily or sub-daily time step. 
%
% Compared to its previous version VarKarst (Hartmann et al., 2015),
% V2Karst:
% (1) uses new definitions for the parameters Vsoi, Vepi and Kepi (they are
%     defined as cell average properties)
% (2) separates the evapotranspiration flux into evaporation from canopy 
%     interception, transpiration and soil evaporation
% (3) comprises three soil layers
% (4) uses the Penman Monteith equation to compute potential
%      evapotranspiration.
% (5) allows for both daily and sub-daily simulations (only daily
%     simulations were possible for VarKarst.
%
% USAGE:
% [Q_epi_d,ET_act_d,Q_surf_d,FLUXES_d,STATES_d,Cont_area_d,time_sim] =...
%               V2Karst(param,n,zm,zh,z,P,Rn,Ta,RH,G,WS,...
%                       time,States_ini,Conc_flow,warmup,Kt)
%
% INPUTS:
%
% PARAMETERS:
%        param = value of model parameters                     - vector(1,15) 
%               param(1) = h_veg [m] (Vegetation height)  
%               param(2) = rs_st [s/m] (Stomatal resistance)
%               param(3) = LAI_min [%](Reduction in leaf area 
%                          index in the dormant season  
%                          compared to the growing season)
%               param(4) = LAI_max [m2/m2] (Annual maximum 
%                          leaf area index during growing
%                          season)
%               param(5) = Vr [mm] (Maximum storage capacity 
%                          of the rooting zone) 
%               param(6) = Vcan [mm per unit of LAI]
%                         (Canopy storage capacity)   
%               param(7) = f_red [-] (reduction factor for 
%                          transpiration in soil layer 3)  
%               param(8) = z0 [m] (soil roughness length)   
%               param(9) = rs_soi [s/m](soil surface 
%                          resistance)
%              param(10) = k [-] (Beer-Lambertâ€™s law 
%                          extinction coefficient)    
%              param(11) = Ve [mm] (Maximum storage 
%                          capacity of first soil layer)
%              param(12) = a [-] (Spatial variability 
%                          parameter)
%              param(13) = Vsoi [mm] (Mean soil storage 
%                        capacity)
%              param(14) = Vepi [mm] (Mean epikarst storage 
%                          capacity)
%              param(15) = Kepi [d] (Mean epikarst outflow 
%                          coefficient)
%            n = number of model vertical compartments [-]         - scalar
%           zm = height of wind speed measurements [m]             - scalar
%           zh = height of humidity measurements [m]               - scalar
%            z = elevation [m]                                     - scalar
% 
% INPUTS DATA: (H is the length of the simulation period)
%           P = daily precipitation                           - vector(H,1)
%          Rn = net radiation [MJ/m2/step]                    - vector(H,1)
%          Ta = mean air temperature                          - vector(H,1)
%               or matrix of mean (Ta(:,1)),                 or vector(H,3)    
%               minimum(Ta(:,2))and maximum(Ta(:,3))
%               temperature [Â°C]   
%              (see Note below)
%         RH = mean relative humidity                         - vector(H,1)
%              or matrix of mean (RH(:,1)),                  or vector(H,3) 
%              minimum (RH(:,2))and maximum(RH(:,3)) 
%              relative humidity [%]   
%             (see Note below)
%          G = ground heat flux [MJ/m2/step]                  - vector(H,1)
%         WS = daily wind speed [m/s]                         - vector(H,1)
%       time = time vector of date numbers (dates at          - vector(H,1)
%              the beginning of the time steps) 
%
% INTIAL STATES:
% States_ini = vector of inital states                        - vector(1,4)
%             (% saturation of deeper vertical compartment)   
%             (see note below)
%              States_ini(1): soil layer 1
%              States_ini(2): soil layer 2
%              States_ini(3): soil layer 3
%              States_ini(4): epikarst layer
%
% OTHER INPUTS:
%  Conc_flow = flag for activation of the concentration flow       - scalar       
%              component of the model                            
%                 - Conc_flow = 1: compartment saturation excess
%                   is channelled to the next unsaturated 
%                   compartment(concentration flow)
%                 - Conc_run = 0: compartment saturation excess
%                   is lost as surface runoff (no concentation flow)
%     warmup = warmup period to be discarded in the simulations    - scalar
%              [months]
%         Kt = time conversion factor [s/time step]                - scalar 
%              for daily time step 3600*24 s/d
%
% OUTPUTS:(H_sim is the length of the simulation period 
%          excluding the warmup period)
%    Q_epi_d = average daily recharge over all            - vector(H_sim,1)
%              compartments [mm]
%   ET_act_d = average daily actual                       - vector(H_sim,1)
%              evapotranspiration over all 
%              compartments (sum of evaporation 
%              from interception, transpiration
%              and soil evaporation)[mm]
%   Q_surf_d = average daily surface runoff over all      - vector(H_sim,1)
%              compartments [mm]
%   FLUXES_d = average daily fluxes over all             - vector(H_sim,12)
%              compartments [mm]
%              FLUXES_d(:,1) = precipitation input 
%              FLUXES_d(:,2) = actual evaporation from 
%                              interception 
%              FLUXES_d(:,3) = actual soil evaporation 
%              FLUXES_d(:,4) = actual transpiration in soil 
%                              layer 1
%              FLUXES_d(:,5) = actual transpiration in soil 
%                              layer 2
%              FLUXES_d(:,6) = actual transpiration in soil 
%                              layer 3
%              FLUXES_d(:,7) = lateral flow (flow from saturated 
%                              to unsaturated compartments)
%              FLUXES_d(:,8) = soil saturation excess
%              FLUXES_d(:,9) = throughfall 
%             FLUXES_d(:,10) = potential evaporation 
%                              from interception
%             FLUXES_d(:,11) = potential soil evaporation
%             FLUXES_d(:,12) = potential transpiration
%   STATES_d = average daily state variables over all     - vector(H_sim,6)
%              compartments [% saturation]
%              STATES_d(:,1) = soil water storage in all 
%                              soil layers
%              STATES_d(:,2) = epikarst water storage
%              STATES_d(:,3) = soil water storage in layers 
%                              1+2(rooting zone)
%              STATES_d(:,4) = soil water storage in layer 1
%              STATES_d(:,5) = soil water storage in layer 2
%              STATES_d(:,6) = soil water storage in layer 3
% Cont_area_d = daily Number of compartments in which     - vector(H_sim,1)
%              the soil generates a saturation excess 
%              flow to the epikarst (contributing areas)[-]
% time_sim = time vector of date numbers for the          - vector(H_sim,1)
%            simulation period (excluding the warmup
%            period)
%
% NOTES:
% - Inputs should verify: Ve < Vr < Vsoi_max(n) where Vsoi_max(n) is the 
%   total storage capacity of the deepest soil compartment.
% - It is recommended to assess the Penman Monteith potential 
%   evapotranspiration using minimum and maximum daily temperature and 
%   minimum and maximum daily relative humidity due to the non-linearity  
%   of the relationships (Allen et al., 1998).
% - Ground heat flux can be neglected at daily time scale (Allen et al.,
%   1998; Shuttleworth, 2012).
% - Wind speed, humidity and temperature must be provided above the
%   vegetation.
% - Initial condition: all compartments initially store the same depth of 
%   water, except for shallow compartments for which the initial condition 
%   is equal to their storage capacity when the initial condition exceed
%   the storage capacity.
% - the seasonality of leaf area index (function LAI_seasonality) is 
%   calculated assuming  LAI is equal to its maximum value in June, July 
%   and August, and to its minimum value in December, January and February.
%
% REFERENCES:
%
% Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998),Crop 
% evapotranspiration: Guidelines for computing crop requirements, FAO 
% Irrigation and Drainage Paper 56, Food and Agriculture Organization 
% (FAO), Rome, Italy.
%
% Bohn, T. J., and E. R. Vivoni (2016), Process-based characterization of 
% evapotranspiration sources over the North American monsoon region, Water 
% Resour. Res., 52(1), 358-384, doi:10.1002/2015WR017934.
%
% Hartmann, A., T. Gleeson, R. Rosolem, F. Pianosi, Y. Wada, and T. Wagener 
% (2015), A large-scale simulation model to assess karstic groundwater 
% recharge over Europe and the Mediterranean, Geosci. Model Dev., 8(6), 
% 1729-1746, doi:10.5194/gmd-8-1729-2015.
% 
% Neitsch, S.L., Arnold, J.G., Kiniry, J.R., Williams, J.R., (2009), 
% Soil and Water Assessment Tool Theoretical Documentation - Version 2009,
% Technical Report no. 406. Texas Water Resources Institute, College  
% Station. Texas. 
%
% Ruiz, L., M. R. R. Varma, M. S. M. Kumar, M. Sekhar, J.-C. Marechal, M. 
% Descloitres, J. Riotte, S. Kumar, C. Kumar, and J.-J. Braun (2010), 
% Water balance modelling in a tropical watershed under deciduous forest 
% (Mule Hole, India): Regolith matric storage buffers the groundwater 
% recharge process, J. Hydrol., 380(3-4), 460-472, 
% doi:10.1016/j.jhydrol.2009.11.020.
%
% Van Dijk, A. I. J. M., and L. A. Bruijnzeel (2001), Modelling rainfall 
% interception by vegetation of variable density using an adapted  
% analytical model. Part 1. Model description, J. Hydrol., 247(3), 230-238, 
% doi:10.1016/S0022-1694(01)00392-4.
%
% This function is part of the V2Karst model V1.1 by F. Sarrazin, A. 
% Hartmann, F. Pianosi, R. Rosolem, T. Wagener (2019, Geosci. Model Dev.)
% V2Karst is provided under the terms of the GNU General Public License 
% version 3.0.
% This function was prepared by Fanny Sarrazin, University of Bristol,
% November 2018 (fanny.sarrazin@bristol.ac.uk).

%--------------------------------------------------------------------------
% 1. Check and recover inputs
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 1.1 Recover parameters
%--------------------------------------------------------------------------
if length(param)~=15 || ~isvector(param)
    error('''param'' must be a vector with length 15')
end
if any(param<0);error('elements in ''param'' must be positive');end
h_veg   = param(1);
rs_st   = param(2);
LAI_min = param(3);
LAI_max = param(4);
Vr      = param(5);
Vcan    = param(6);
f_red   = param(7);
z0      = param(8);
rs_soi  = param(9);
k       = param(10);
Ve      = param(11);
a       = param(12);
Vsoi    = param(13);
Vepi    = param(14);
Kepi    = param(15);
if ~isscalar(n) || n<0; error('''n'' must be scalar and positive');end
if mod(n,1)~=0;error('''n'' must be a integer');end
if ~isscalar(z) || z<0; error('''z'' must be scalar and positive');end
if ~isscalar(zm); error('''zm'' must be scalar and positive');end
if ~isscalar(zh); error('''zh'' must be scalar and positive');end
if h_veg>zm || h_veg>zh;error('h_veg (param(1)) must be lower than zm and zh');end

%--------------------------------------------------------------------------
% 1.2 Check forcing data
%--------------------------------------------------------------------------
if size(time,2)~=1 || size(P,2)~=1 || size(Rn,2)~=1 || size(G,2)~=1 ||...
        size(WS,2)~=1 
    error('''time'',''P'',''Rn'',''G'' and''WS'' must be column vectors')
end
if ~any(size(Ta,2)==[1 3]) || ~any(size(RH,2)==[1 3])
    error('''Ta'' and ''RH'' must have 1 or 3 columns')
end
H=size(P,1); % length of simulation period
if size(time,1)~=H || size(Rn,1)~=H || size(G,1)~=H || size(WS,1)~=H ||...
        size(Ta,1)~=H ||size(RH,1)~=H 
    error('''time'',''P'',''Rn'',''G'',''WS'',''Ta''and ''RH'' must have the same number of rows')
end

%--------------------------------------------------------------------------
% 1.3 Check and recover initial states
%--------------------------------------------------------------------------
if length(States_ini)~=4 || ~isvector(States_ini)
    error('''States_ini'' must be a vector with length 4')
end

if any(States_ini>100 | States_ini<0)
    error('elements in ''States_ini'' must be between 0 and 100')
end
Vsoi1_ini = States_ini(1);
Vsoi2_ini = States_ini(2);
Vsoi3_ini = States_ini(3);
Vepi_ini = States_ini(4);

%--------------------------------------------------------------------------
% 1.3 Check other inputs
%--------------------------------------------------------------------------
if ~isscalar(Conc_flow);error('''Conc_flow'' must be scalar');end
if ~any(Conc_flow==[0 1]);error('''Conc_flow'' must be equal to 0 or 1');end
    
if ~isscalar(warmup);error('''warmup'' must be scalar');end
if mod(warmup,1)~=0;error('''warmup'' must be a integer');end
if warmup<0;error('''warmup'' must be a positive');end
if Kt<0 || Kt-round(Kt)~=0;error('Kt should be a positive integer');end
if Kt>86400; error('Kt must be below or equal to 86400 (simulations can be daily or subdaily');end

%--------------------------------------------------------------------------
% 2. Define constant parameters
%--------------------------------------------------------------------------
d = 0; % zero plane displacement height for bare soil [m]

%--------------------------------------------------------------------------
% 3. Seasonality of leaf area index 
%--------------------------------------------------------------------------
 LAI = LAI_seasonality(LAI_max,LAI_min,time); % daily cell average leaf 
                                              % area index [m2/m2]

%--------------------------------------------------------------------------
% 4. Vegetation cover fraction
%--------------------------------------------------------------------------
% Vegetation cover fraction is assessed using the Beer-Lambert's law as in
% [van Dijk and Bruijnzeel, 2001; Ruiz et al., 2010]
fc = 1-exp(-k*LAI); % daily vegetation cover fraction [-]

%--------------------------------------------------------------------------
% 5. Calculate potential evapotranspiration rates
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 5.1 Atmospheric pressure (Allen et al., 1998, Eq.(7))
%--------------------------------------------------------------------------
Pa = 101.3*((293-0.0065*z)/293)^5.26; %[kPa]

%--------------------------------------------------------------------------
% 5.2 Canopy and soil aerodynamic resistances
%--------------------------------------------------------------------------
% Canopy roughness lengths used in Sarrazin et al. (2017)
%(Allen et al., 1998, BOX 4, p.21)
zom_can = 0.123*h_veg; % canopy roughness length for momentum transfer [m]
zoh_can = 0.1*zom_can; % canopy roughness length for heat and vapor 
                       % transfer [m] 

% Other option to assess canopy roughness lengths
% (Neitsch et al., 2009, p128)
% if h_veg > 2
%     zom_can=(0.058*(h_veg*100)^1.19)/100; 
% else
%     zom_can = 0.123*h_veg;
% end
% zoh_can = 0.1*zom_can;

% Canopy zero plane displacement height(Allen et al., 1998, BOX 4, p.21)
d_can = 2/3*h_veg; %  [m]  

% Canopy aerodynamic resistance [s/m]
ra_can = aerodynamic_resistance(zom_can,zoh_can,d_can,WS,zm,zh);
% Bare soil aerodynamic resistance [s/m]
ra_soi = aerodynamic_resistance(z0,z0,d,WS,zm,zh);

%--------------------------------------------------------------------------
% 5.3 Canopy surface resistance
%--------------------------------------------------------------------------
% Effective canopy storage capacity (cell average LAI is rescaled to
% the vegetated fraction following Bohn and Vivoni, 2016)
rs_can = rs_st./(LAI./fc);

%--------------------------------------------------------------------------
% 5.4 Potential evapotranspiration rates
%--------------------------------------------------------------------------
% Potential canopy interception [mm]
Ecan_pot = Penman_Monteith(0,ra_can,Rn,Ta,RH,G,Pa,Kt); 
% Potential transpiration [mm]
T_pot = Penman_Monteith(rs_can,ra_can,Rn,Ta,RH,G,Pa,Kt); 
% Potential soil evaporation [mm]
Es_pot = Penman_Monteith(rs_soi,ra_soi,Rn,Ta,RH,G,Pa,Kt);

%--------------------------------------------------------------------------
% 6. Evaporation from canopy interception
%--------------------------------------------------------------------------
[Tf,t_wet,Ecan_act] = interception_routine(Vcan,LAI,fc,P,Ecan_pot,Kt);

%--------------------------------------------------------------------------
% 7. Soil and epikarst water balance
%--------------------------------------------------------------------------
[Q_epi_avg,ETsoi_act_avg,Q_surf_avg,STATES,FLUXES,Cont_area] =...
    soil_epikarst_routine(fc,f_red,Vr,Ve,a,Vsoi,Vepi,Kepi,n,Tf,T_pot,...
    Es_pot,t_wet,Vsoi1_ini,Vsoi2_ini,Vsoi3_ini,Vepi_ini,Conc_flow);

%--------------------------------------------------------------------------
% 8. Remove warmup period
%--------------------------------------------------------------------------
% Prepare time vectors
date_vec = datevec(time);
Year_month = unique(date_vec(:,1:2),'rows'); % months of the simulation 
                                             % period (including warmup)
Year_month_warm = Year_month(warmup+1:end,:); % months of the simulation 
                                              % period (excluding warmup)
% First day following the warmup period  
idx_start = find(ismember(date_vec(:,1:4),[Year_month_warm(1,:) 1 0],'rows')>0);
% Time vector for simulation period excluding warmup
time_sim=time(idx_start:end);

% Assess daily outputs
Q_epi_d  = Q_epi_avg(idx_start:end);
ET_act_d = Ecan_act(idx_start:end)+ETsoi_act_avg(idx_start:end);
Q_surf_d = Q_surf_avg(idx_start:end);
FLUXES_d = [P(idx_start:end),Ecan_act(idx_start:end),...
           FLUXES(idx_start:end,:),Tf(idx_start:end),...
           Ecan_pot(idx_start:end),Es_pot(idx_start:end),...
           T_pot(idx_start:end)];
STATES_d = STATES(idx_start:end,:);
Cont_area_d = Cont_area(idx_start:end);


