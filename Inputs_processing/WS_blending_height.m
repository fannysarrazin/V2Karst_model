function WS_b = WS_blending_height(zom,d,WS,zm,zb)

% This function computes wind speed at a blending height following 
% Lhomme et al.(2014).
%
% USAGE:
%  WS_b = WS_blending_height(zom,d,WS,zm,zb)
%
% INPUTS:
% zom = roughness length for momentum transfer [m]             - scalar
%   d = zero plane displacement height [m]                     - scalar              
%  WS = wind speed measurement[m/s]                            - vector(H,1)
%  zm = height of wind speed measurements [m]                  - scalar
%  zb = blending heigth [m]                                    - scalar
%
% OUTPUTS:
% WS_b = wind speed at blending height [m s-1]                - vector(H,1)
%
% REFERENCE:
% Lhomme, J. P., N. Boudhina, and M. M. Masmoudi (2014), Technical Note: 
% On the Matt–Shuttleworth approach to estimate crop water requirements, 
% Hydrol. Earth Syst. Sci., 18(11), 4341–4348, doi:10.5194/hess-18-4341-2014.

% Calculate wind speed at blending height (Lhomme et al., 2014, Eq.(6))
WS_b = WS*log((zb-d)/zom)./log((zm-d)/zom);

% Check variable
if any(isnan(WS_b)); error('WS_b contains NaN');end
if any(WS_b<0);error('WS_b<0');end
if any(WS_b>300);error('WS_b>300');end
    