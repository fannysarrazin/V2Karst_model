function [Q_epi_T,ET_act_T,Q_surf_T,FLUXES_T,STATES_T,...
          Cont_area_T,time_sim,Delta_V_T] =...
          V2Karst(param,n,zm,zh,Pa,P,Rad,Ta,RH,G,WS,Swe,CO2,time,...
          States_ini,Conc_flow,warmup,Kt)

% This function simulates the V2Karst model, which is a vegetation-
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
%      evapotranspiration
% (5) allows for both daily and sub-daily simulations (only daily
%     simulations were possible for VarKarst).
%
% USAGE:
% [Q_epi_T,ET_act_T,Q_surf_T,FLUXES_T,STATES_T,Cont_area_T,time_sim,
%  Delta_V_T] = V2Karst(param,n,zm,zh,Pa,P,Rad,Ta,RH,G,WS,Swe,CO2,time,...
%                       States_ini,Conc_flow,warmup,Kt)
%h_veg     = param(1);
% 
% we denote T the simulation time step.
% 
% INPUTS:
%
% PARAMETERS:
%        param = value of model parameters                   - vector(1,15) 
%               param(1) = h_veg [m] (Vegetation height)  
%               param(2) = rs_st [s m-1] (Stomatal resistance)
%               param(3) = LAI_min [%] (Reduction in leaf area 
%                          index in the dormant season  
%                          compared to the growing season)
%                          (see Note 1)
%               param(4) = LAI_max [m2 m-2] (Annual maximum 
%                          leaf area index during the growing
%                          season)
%               param(5) = Vr [mm] (Maximum storage capacity 
%                          of the rooting zone) 
%                          (see Note 2 below)
%               param(6) = Vcan [mm per unit of LAI]
%                         (Canopy storage capacity)   
%               param(7) = k [-] (Beer-Lambert’s law 
%                          extinction coefficient) 
%               param(8) = f_red [-] (reduction factor for 
%                          transpiration in soil layer 3)  
%               param(9) = z0 [m] (soil roughness length)   
%              param(10) = rs_soi [s m-1](soil surface resistance)     
%              param(11) = Ve [mm] (Maximum storage 
%                          capacity of first soil layer)
%                          (see Note 2 below)
%              param(12) = a [-] (Spatial variability 
%                          parameter)
%              param(13) = Vsoi [mm] (Mean soil storage 
%                          capacity)
%                          (see Note 2 below)
%              param(14) = Vepi [mm] (Mean epikarst storage 
%                          capacity)
%              param(15) = Kepi [T-1] (Mean epikarst outflow 
%                          coefficient)
%              param(16) = S_rs [ppm-1] (sensitivity of 
%                         stomatal  resistance to CO2)   
%                         (by setting S_rs=0 the sensitivity
%                         of stomatal resistance ot CO2 is 
%                         is removed
%
% Optional parameters (see Note 3 below): 
%            param(17) = alpha_veg [-] (vegetation albedo) 
%            param(18) = alpha_soi [-] (bare soil albedo) 
%
%            n = number of model vertical compartments [-]         - scalar
%           zm = height of wind speed measurements [m]             - scalar
%           zh = height of humidity measurements [m]               - scalar
% 
% INPUTS DATA: (H is the length of the simulation period)
%          Pa = atmospheric pressure[kPa]                     - vector(H,1)
%               (see Note 4 below)
%           P = precipitation [mm T-1]                        - vector(H,1)
%         Rad = radiation input [MJ m-2 T-1]                  - vector(H,1)
%               Net radiation is either provided as input  or - vector(H,2)
%               or calculated based on downward radiation.
%               - option 1: vector of size(H,1) that       
%                 contains values of net radiation
%               - option 2: vector of size(H,2) that contains 
%                 values of downward radiation. With this 
%                 option, net radiation will be computed in
%                 V2Karst based on albedo and air temperature
%                 as a proxy for surface temperature 
%                 (see help of net_radiation.m for further 
%                  explanation).
%                    *Rad(:,1): short wave downward radiation
%                    *Rad(:,2): long wave downward radiation
%               (see Note 5 below)
%          Ta = mean air temperature [deg C]                  - vector(H,1)
%               or matrix of mean (Ta(:,1)),                 or vector(H,3)    
%               minimum(Ta(:,2))and maximum(Ta(:,3)) [deg C] 
%              (see Note 6 and 7 below)
%         RH = mean relative humidity [%]                     - vector(H,1)
%              or matrix of mean (RH(:,1)),                  or vector(H,3) 
%              minimum (RH(:,2))and maximum(RH(:,3)) [%]
%             (see Note 6 and 7 below)
%          G = ground heat flux [MJ m-2 T-1]                  - vector(H,1)
%              (see Note 8 below)
%         WS = daily wind speed [m s-1]                       - vector(H,1)
%              (see Note 7 below)   
%        Swe = daily snow water equivalent [mm]               - vector(H,1)
%        CO2 = atmospheric CO2 concentration [ppm]            - vector(H,1)
%              CO2 is either a time series or a scalar        or scalar
%              in this case it is assumed to be constant 
%              over time)
%       time = time vector of date numbers (dates at          - vector(H,1)
%              the beginning of the time steps) 
%
% INTIAL STATES:
% States_ini = vector of inital states [%]                    - vector(1,4)
%             (% saturation of deeper vertical compartment)   
%             (see Note 9 below)
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
%         Kt = time step in seconds [s T-1]                        - scalar 
%              for daily time step set to 3600*24 s d-1
%
% OUTPUTS:(H_sim is the length of the simulation period 
%          excluding the warmup period)
%    Q_epi_T = time series of average recharge over all   - vector(H_sim,1)
%              compartments [mm T-1]
%   ET_act_T = time series of average actual              - vector(H_sim,1)
%              evapotranspiration over all 
%              compartments (sum of evaporation 
%              from interception, transpiration
%              and soil evaporation)[mm T-1]
%   Q_surf_T = time series of average surface runoff     - vector(H_sim,1)
%              over all compartments [mm T-1]
%   FLUXES_T = time series of average fluxes over all    - vector(H_sim,12)
%              compartments [mm T-1]
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
%   STATES_d = time series of average state variables     - vector(H_sim,6)
%              over all  compartments [% saturation]
%              STATES_d(:,1) = soil water storage in all 
%                              soil layers
%              STATES_d(:,2) = epikarst water storage
%              STATES_d(:,3) = soil water storage in layers 
%                              1+2(rooting zone)
%              STATES_d(:,4) = soil water storage in layer 1
%              STATES_d(:,5) = soil water storage in layer 2
%              STATES_d(:,6) = soil water storage in layer 3
% Cont_area_d = time series of number of compartments in  - vector(H_sim,1)
%              which the soil generates a saturation excess 
%              flow to the epikarst (contributing areas)[-]
%    time_sim = time vector of date numbers for the          - vector(H_sim,1)
%               simulation period (excluding the warmup
%               period)
%
% NOTES:
% 1. The seasonality of leaf area index (function LAI_seasonality) is 
%    calculated assuming  LAI is equal to its maximum value in June, July 
%    and August, and to its minimum value in December, January and
%    February.
% 2. Inputs should verify: Ve < Vr < Vsoi_max(n) where Vsoi_max(n) is the 
%    total storage capacity of the deepest soil compartment.
% 3. alpha_veg and alpha_soi need to be specified when net radiation is not 
%    directly provided as input but calculated from downward longwave and 
%    shortwave radiations. 
% 4. When Pa is not available, it can be estimated based on elevation as 
%    reported in Allen et al., 1998, Eq.(7) (this calculation has to be 
%    performed prior to calling V2Karst).
% 5. Net radiation in assessed from downward radiations assuming that air 
%    temperature is a proxy for surface temperature, which has only been 
%    tested for daily time step. When net radiation has to be calculated
%    from downward radiations, we recommend to use a daily time step (and
%    NOT sub-daily) to simulate V2Karst.
% 6. For daily time step, it is recommended to assess the Penman
%    Monteith potential evapotranspiration using minimum and maximum daily 
%    temperature and minimum and maximum daily relative humidity due to the 
%    non-linearity of the relationships (Allen et al., 1998).
% 7. Wind speed, relative humidity and temperature must be provided above 
%    the vegetation. If this is not the case, relative humidity and wind
%    speed can be calculated at a blending height using the function
%    Inputs_processing/RH_WS_blending_height.m
% 8. Ground heat flux can be neglected at daily time scale (Allen et al.,
%    1998; Shuttleworth, 2012).
% 9. Initial condition: all compartments initially store the same depth of 
%    water, except for shallow compartments for which the initial condition 
%    is equal to their storage capacity when the initial condition exceed
%    the storage capacity.

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
% Yang, Y., M. L. Roderick, S. Zhang, T. R. McVicar and R. J. Donohue
% (2019), Hydrologic implications of vegetation response to elevated CO2 in 
% climate projections, Nature Climate Change, 9, 44-48, 
% doi:10.1038/s41558-018-0361-0
%
% Van Dijk, A. I. J. M., and L. A. Bruijnzeel (2001), Modelling rainfall 
% interception by vegetation of variable density using an adapted  
% analytical model. Part 1. Model description, J. Hydrol., 247(3), 230-238, 
% doi:10.1016/S0022-1694(01)00392-4.
%
% This function is part of the V2Karst model by F. Sarrazin, A. Hartmann, 
% F. Pianosi, R. Rosolem, T. Wagener (2018, Geosci. Model Dev.)
% V2Karst is provided under the terms of the GNU General Public License 
% version 3.0.
% This function was prepared by Fanny Sarrazin (fanny.sarrazin@ufz.de).

%--------------------------------------------------------------------------
% 1. Check and recover inputs
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 1.1 Recover parameters
%--------------------------------------------------------------------------
if ~isvector(param)
    error('''param'' must be a vector')
end

% Set default value of optional parameters alpha_veg and alpha_soi: 
alpha_veg = [];
alpha_soi = [];

if any(param<0);error('elements in ''param'' must be positive');end
if length(param) == 18
    alpha_veg = param(17);
    alpha_soi = param(18);
elseif length(param)>18 || length(param)<16
    error('''param'' must be a vector with length between 15 and 18')
end

h_veg     = param(1);
rs_st     = param(2);
LAI_min   = param(3);
LAI_max   = param(4);
Vr        = param(5);
Vcan      = param(6);
k         = param(7);
f_red     = param(8);
z0        = param(9);
rs_soi    = param(10);
Ve        = param(11);
a         = param(12);
Vsoi      = param(13);
Vepi      = param(14);
Kepi      = param(15);
S_rs      = param(16);

if ~isscalar(n) || n<0; error('''n'' must be scalar and positive');end
if mod(n,1)~=0;error('''n'' must be a integer');end
if ~isscalar(zm); error('''zm'' must be scalar and positive');end
if ~isscalar(zh); error('''zh'' must be scalar and positive');end
if h_veg>zm || h_veg>zh;error('h_veg (param(1)) must be lower than zm and zh');end

%--------------------------------------------------------------------------
% 1.2 Check forcing data
%--------------------------------------------------------------------------
if size(time,2)~=1 || size(P,2)~=1 || size(G,2)~=1 ||...
        size(WS,2)~=1 ||size(Pa,2)~=1 ||size(Swe,2)~=1   
    error('''time'',''P'',''Rn'',''G'',''WS'',''Pa'' and ''Swe'' must be column vectors')
end
if ~any(size(Ta,2)==[1 3]) || ~any(size(RH,2)==[1 3])
    error('''Ta'' and ''RH'' must have 1 or 3 columns')
end
H=size(P,1); % length of simulation period
if size(time,1)~=H || size(Rad,1)~=H || size(G,1)~=H || size(WS,1)~=H ||...
        size(Ta,1)~=H ||size(RH,1)~=H || size(Swe,1)~=H ||size(Pa,1)~=H
    error('''time'',''P'',''Rn'',''G'',''WS'',''Ta'',''RH'',''Pa'' and ''Swe'' must have the same number of rows')
end

if size(Rad,2)>2
    error('''Rad'' must have 1 or 2 columns')
elseif size(Rad,2)==2 && isempty(alpha_veg)
    error('''alpha_veg'' and ''alpha_soi'' need to be specified to compute net radiation (''param'' must have 18 elements)')
end

% Check CO2 input:
if ~isscalar(CO2) && length(CO2)~=H
    error('CO2 must be scalar or it must be a column vector with the same number of rows as P')
end
if size(CO2,2)~=1
    error('''CO2'' be column vectors')
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
if mod(warmup,1)~=0;error('''warmup'' must be integer');end
if warmup<0;error('''warmup'' must be positive');end
if Kt<0 || Kt-round(Kt)~=0;error('Kt should be a positive integer');end
if Kt>86400; error('Kt must be below or equal to 86400 (simulations can be daily or subdaily');end

%--------------------------------------------------------------------------
% 2. Define constant parameters
%--------------------------------------------------------------------------
d = 0; % zero plane displacement height for bare soil [m]
emissivity = 1; % [-] surface emissivity (to calculate net radiation from 
% downward radiation)
swe_thres = 0.1; % [mm] threshold on swe above which we consider the presence
% of the snow pack and set the potential transpiration and soil evaporation
% to 0.

%--------------------------------------------------------------------------
% 3. Seasonality of leaf area index 
%--------------------------------------------------------------------------
 LAI = LAI_seasonality(LAI_max,LAI_min,time); % cell average leaf area index 
                                              % [m2/m2]

%--------------------------------------------------------------------------
% 4. Vegetation cover fraction
%--------------------------------------------------------------------------
% Vegetation cover fraction is assessed using the Beer-Lambert's law as in
% [van Dijk and Bruijnzeel, 2001; Ruiz et al., 2010]
fc = 1-exp(-k*LAI); % vegetation cover fraction [-]

%--------------------------------------------------------------------------
% 5. Calculate potential evapotranspiration rates
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 5.1 Net radiation 
%--------------------------------------------------------------------------
if size(Rad,2)==2
    % Net radiation is assessed from downward radiation using the average 
    % albedo over the vegetated and non-vegetated fracions. Use mean temperature
    Rn = net_radiation(alpha_veg*fc+alpha_soi*(1-fc),emissivity,Rad(:,1),Rad(:,2),Ta(:,1));
else
    % Net radiation is directly provided as input
    Rn = Rad;
end

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
rs_can = rs_st./(LAI./fc).*(1+S_rs.*(CO2-300));
% The term (1+S_rs*(CO2-300)) captures the variability of the surface
% resistance depending on the atmospheric CO2 concentration (Yang et al.,
% 2019)

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
% 7. Snow routine
%--------------------------------------------------------------------------
% Assess effective precipitation (Pre_eff) from the variation of the
% snow pack and set to 0 potential transpiration and soil evaporation in
% presence of a snow pack.
[Pre_eff,T_pot_eff,Es_pot_eff,Swe_new] =...
                              snow_routine(Swe,Tf,T_pot,Es_pot,swe_thres);

%--------------------------------------------------------------------------
% 8. Soil and epikarst water balance
%--------------------------------------------------------------------------
[Q_epi_avg,ETsoi_act_avg,Q_surf_avg,STATES,FLUXES,Cont_area,Delta_V] =...
    soil_epikarst_routine(fc,f_red,Vr,Ve,a,Vsoi,Vepi,Kepi,n,Pre_eff,T_pot_eff,...
    Es_pot_eff,t_wet,Vsoi1_ini,Vsoi2_ini,Vsoi3_ini,Vepi_ini,Conc_flow);

%--------------------------------------------------------------------------
% 9. Remove warmup period
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

% Output time series
Q_epi_T  = Q_epi_avg(idx_start:end);
ET_act_T = Ecan_act(idx_start:end)+ETsoi_act_avg(idx_start:end);
Q_surf_T = Q_surf_avg(idx_start:end);
FLUXES_T = [P(idx_start:end),Ecan_act(idx_start:end),...
           FLUXES(idx_start:end,:),Tf(idx_start:end)];
STATES_T = STATES(idx_start:end,:);
Cont_area_T = Cont_area(idx_start:end);
Delta_V_T = [Delta_V(end)+Swe_new(end)-Swe_new(1),Delta_V(end)-Delta_V(idx_start-1)+Swe_new(end)-Swe_new(idx_start-1)];


%--------------------------------------------------------------------------
% 10. Check water balance
%--------------------------------------------------------------------------
if abs(sum(P - Q_epi_avg - Ecan_act - ETsoi_act_avg - Q_surf_avg) - Delta_V_T(1))>10^(-5)
    error('Water balance not satisfied')
end


    
    