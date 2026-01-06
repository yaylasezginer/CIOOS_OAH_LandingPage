function [fCO2,pH,CO2,HCO3,CO3,omega_cal,omega_arag] = eh_carbeq(S,T,P,DIC,Alk)

% eh_carbeq   Concentrations of carbonate species
%=========================================================================
% eh_carbeq  Version 1.5 (June 14, 2010)
%
% USAGE:  [fCO2,pH,CO2,HCO3,CO3,omega_cal,omega_arg] = eh_carbeq(S,T,P,DIC,Alk)
%
% DESCRIPTION:
%    Concentrations of carbonate species including fCO2 and pH from DIC and
%       Alkalinity as a function of temperature, salinity, and pressure.
%    Also calculates saturation state for calcite and aragonite.
%    Considers only carbon, boron, and water ions in alkalinity, using 
%       equations from Appendix 5A.1 (b) of Emerson and Hamme (2022).
%    This function requires simpler inputs than CO2SYS and is suitable for
%       many educational applications. A comparison between results of this
%    program and CO2SYS can be made by calculating the carbonate species for
%    the conditions in Table 5.5 of E&H(22).  For surface waters they are
%    identical, but for deeper waters they are slightly different because of 
%   the silicic acid and phosphate contributions
%
% INPUT:  (if inputs are not singular they must have same dimensions)
%   S = practical salinity    [PSS-78 scale]
%   T = temperature [degree C, ITS-90 scale]
%   P = pressure [dbar, use 0 for surface]
%   DIC = Dissolved Inorganic Carbon   [umol/kg]
%   Alk = Alkalinity   [ueq/kg]
%
% OUTPUT:
%   fCO2 = fugacity of CO2  [uatm]
%   pH = -log10 of hydrogen ion concentration [total pH scale]
%   CO2 = concentration of dissolved CO2 and H2CO3  [micromol/kg]
%   HCO3 = concentration of dissolved HCO3  [micromol/kg]
%   CO3 = concentration of dissolved CO3  [micromol/kg]
%   omega_cal = degree of saturation for calcite (omega < 1 is
%   undersaturated)
%   omega_arag = degree of saturation for aragonite (omega < 1 is
%   undersaturated)
% 
% AUTHOR:  Roberta Hamme (University of Victoria) rhamme@uvic.ca
%
% REFERENCE:
%    This program uses equilibrium constants and equations from:
%    Dickson, A.G., Sabine, C.L. and Christian, J.R. (Eds.) 2007. Guide to
%    best practices for ocean CO2 measurements. PICES Special Publication 
%    3, 191 pp.
%    Calcite and aragonite solubility from Mucci "The solubility of calcite 
%    and aragonite in seawater at various salinities, temperatures, and one 
%    atmosphere total pressure" 1983. Am. J. Sci., 283. 
%    Mimization routine and pressure corrections from Lewis, E., D.W.R 
%    Wallace, CO2SYS: Program developed for CO2 system calculations
%
% VERSION 1.5 : 25 March 2022
% AUTHOR: Roberta C. Hamme (University of Victoria) 
% This software is available from http://www.cambridge.org/emerson-hamme
% as part of Chemical Oceanography: Element Fluxes in the Sea (2022) 
% by Steven R. Emerson and Roberta C. Hamme
%=========================================================================

%----------------------
% Check input parameters
%----------------------

% check number of input parameters is correct
if nargin ~=5
   error('Must pass 5 input parameters')
end %if

% check datatypes of input parameters are correct
validateattributes(S,{'numeric'},{'nonempty'},mfilename,'S:salinity')
validateattributes(T,{'numeric'},{'nonempty'},mfilename,'T:temperature')
validateattributes(P,{'numeric'},{'nonempty'},mfilename,'P:pressure')
validateattributes(DIC,{'numeric'},{'nonempty'},mfilename,'DIC:Dissolved Inorganic Carbon')
validateattributes(Alk,{'numeric'},{'nonempty'},mfilename,'Alk:Alkalinity')

% Check inputs dimension and verify they have the same shape or are singular
[rs,cs] = size(S);
[rt,ct] = size(T);
[rp,cp] = size(P);
[rd,cd] = size(DIC);
[ra,ca] = size(Alk);

if (rs+cs)>2
    rall=rs; call=cs;
end
if (rt+ct)>2
    if exist('rall','var')
        if (rt ~= rall || ct ~= call)
            error('inputs must have same dimensions or be singular')
        end
    else rall=rt; call=ct;
    end
end
if (rp+cp)>2
    if exist('rall','var')
        if (rp ~= rall || cp ~= call)
            error('inputs must have same dimensions or be singular')
        end
    else rall=rp; call=cp;
    end
end
if (rd+cd)>2
    if exist('rall','var')
        if (rd ~= rall || cd ~= call)
            error('inputs must have same dimensions or be singular')
        end
    else rall=rd; call=cd;
    end
end
if (ra+ca)>2 && exist('rall','var') && (ra ~= rall || ca ~= call)
    error('inputs must have same dimensions or be singular')
end

%----------------------
% Do some basic unit conversion
%----------------------

% calculate temperature in Kelvin
TK = T+273.15;
% gas constant (cm^3 bar^-1 mol^-1 K^-1)
R = 83.1446;
% convert pressure from dbar to bar
P = 0.1*P;
% convert DIC to mol/kg
DIC = DIC * .000001;
% Convert Alkalinity to eq/kg
Alk = Alk * .000001;    

%----------------------
% Calculate ion concentrations and equilibrium constants
%----------------------

% Calculate total borate (BT) from salinity 
% (Uppström, L., Deep-Sea Research, 21,161-162, 1974)
BT = 0.0004157 .* S ./ 35;        % [mol kg^-1]

% Calculate total calcium from salinity 
% (Riley, R.F., and M. Tongudai, Chem. Geol., 2, 263–269, 1967)
Ca = 0.0102846 .* S ./ 35;          % [mol kg^-1]

% Calculate Henry's Law coeff for CO2 (KH) from temp & sal 
% (Weiss, R.F., Mar. Chem., 2, 203-215, 1974)
KH = exp(-60.2409 + 9345.17./TK + 23.3585*log(TK/100) + S.*(0.023517 - 0.00023656*TK + 0.0047036*(TK/100).^2)); % [mol kg^-1 atm^-1]

% Calculate borate equil constant (KB) from temp & sal
% (Dickson, A.G., 1990, Deep-Sea Res., 37, 755-766)
KB_0dbar = exp((-8966.9 - 2890.53*S.^0.5 - 77.942*S + 1.728*S.^1.5 - 0.0996*S.^2)./TK + 148.0248 + 137.1942*S.^0.5 + 1.62142*S - (24.4344 + 25.085*S.^0.5 + 0.2474*S).*log(TK) + 0.053105*S.^0.5.*TK);
KB = KB_0dbar .* exp((-(-29.48 + 0.1622*T - 2.608e-3*T.^2) + 0.5*(-2.84e-3).*P).*P./(R*TK));

% Calculate carbonate equil constants (K1 & K2) from temp & sal 
% (Lueker, T.J., A.G. Dickson, C.D. Keeling, 2000, Mar. Chem., 70, 105-119)
% these are on the total pH scale
K1_0dbar = 10.^(-(3633.86./TK - 61.2172 + 9.6777*log(TK) - 0.011555*S + 0.0001152*S.^2));
K1 = K1_0dbar .* exp((-(-25.5 + 0.1271*T) + 0.5*(-3.08e-3 + 8.77e-5*T).*P).*P./(R*TK));
K2_0dbar = 10.^(-(471.78./TK + 25.929 - 3.16967*log(TK) - 0.01781*S + 0.0001122*S.^2));
K2 = K2_0dbar .* exp((-(-15.82 - 0.0219*T) + 0.5*(1.13e-3 + -1.475e-4*T).*P).*P./(R*TK));

% Calculate water equil constant (KW) from temp & sal 
% (Millero, F.J., 1995, Geochim. Cosmochim. Acta, 59, 661-677)
KW_0dbar = exp(148.9652 - 13847.26./TK - 23.6521*log(TK) + (118.67./TK - 5.977 + 1.0495.*log(TK)).*S.^0.5 - 0.01615*S);
KW = KW_0dbar .* exp((-(-25.60 + 0.2324*T - 3.6246e-3*T.^2) + 0.5*(-5.13e-3 + 7.94e-5*T).*P).*P./(R*TK));

% Calculate Ksp for calcite and argonite 
Ksp_calcite_0dbar = 10.^(-(171.9065 + 0.077993*TK - 2839.319./TK - 71.595*log10(TK) + (0.77712 - 0.0028426*TK - 178.34./TK).*S.^0.5 + 0.07711*S - 0.0041249*S.^1.5));
Ksp_calcite = Ksp_calcite_0dbar .* exp((-(-48.76 + 0.5304*T) + 0.5*(-1.176e-2 + 3.692e-4*T).*P).*P./(R*TK));
Ksp_aragonite_0dbar = 10.^(-(171.945 + 0.077993*TK - 2903.293./TK - 71.595*log10(TK) + (0.068393 - 0.0017276*TK - 88.135./TK).*S.^0.5 + 0.10018*S - 0.0059415*S.^1.5));
Ksp_aragonite = Ksp_aragonite_0dbar .* exp((-(-45.96 + 0.5304*T) + 0.5*(-1.176e-2 + 3.692e-4*T).*P).*P./(R*TK));        

%----------------------
% Solve for H ion concentration using a minimization routine
% based on Newton's method following CO2SYS code
%----------------------

pH = 8;            % initial guess for minimization is pH = 8;
pHTol = 0.0001;     % final pH estimate must be within 0.0001 of true value
deltapH = pHTol+1;  % give an initially large deltapH to start the loop
while any(abs(deltapH) > pHTol)
    H = 10.^(-pH);
    HCO3 = DIC .* K1 .* H ./ (H.^2 + K1 .* H + K1 .* K2);
    CO3 = DIC .* K1 .* K2 ./ (H.^2 + K1 .* H + K1 .* K2);
    BOH4 = KB .* BT ./ (H + KB);
    OH = KW ./ H;
    AlkResid  = Alk - HCO3 - 2*CO3 - BOH4 - OH + H;
    % find Slope dTA/dpH; directly from CO2SYS code
    % (this is not exact, but keeps all important terms);
    Slope = log(10).*(DIC.*K1.*H.*(H.*H + K1.*K2 + 4.*H.*K2)./((H.^2 + K1.*H + K1.*K2).^2) + BOH4.*H./(KB + H) + OH + H);
    deltapH = AlkResid./Slope;
    % to keep the jump from being too big;
    while any(abs(deltapH) > 0.5)
        index = abs(deltapH) > 0.5; 
        deltapH(index) = deltapH(index)./2;
    end %end step checker while
    pH = pH + deltapH;
end %end minimization while

%----------------------
% calculate/convert units HCO3, CO3 and CO2, fCO2 and pH
%----------------------

HCO3 = HCO3 * 1e6;
CO3 = CO3 * 1e6;
CO2 = DIC./(1 + K1./H + K1.*K2./(H.*H)) * 1e6;
fCO2 = CO2 ./ KH;
omega_cal = Ca .* (CO3/1e6) ./ Ksp_calcite;
omega_arag = Ca .* (CO3/1e6) ./ Ksp_aragonite;