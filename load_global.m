% ***********************************************************************************
%                          L O A D  G L O B A L  for
%            C A R D I O V A S C U L A R   S Y S T E M S   M O D E L
% ***********************************************************************************
% Copyright (C) 2019 N. P. Woodall, K. G. Kim, A. L. Colunga, T. F. Dardas, 
% M. S. Olufsen, and  B. E. Carlson
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software"), to deal in
% the Software without restriction, including without limitation the rights to
% use, copy, modify, and merge the Software subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% The authors be cited in any work resulting from the use of the Software.
% The associated published article arXiv:1812.11857
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, OR DATA, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
% ***********************************************************************************
% Set nominal parameters, initial conditions, and upper and lower bounds of ODE
% model
% ***********************************************************************************
function [pars, Init, hi, low, gPars] = load_global(data)

gPars.ABS_TOL  = 1e-12;
gPars.REL_TOL  = gPars.ABS_TOL/10;
gPars.DIFF_INC = sqrt(gPars.ABS_TOL);

% Right Ventricle 
P_RVsyst   = data.P_RVsyst;                       % Systolic RV pressure (mmHg)     *
P_RVdiast  = data.P_RVdiast;                      % Diastolic RV pressure (mmHg)    *

% Pulm Art
P_PAsyst   = data.P_PAsyst;                       % Systolic pulm art press (mmHg)  *
P_PAdiast  = data.P_PAdiast;                      % Diast pulm art press (mmHg)     *
P_PApp     = P_PAsyst - P_PAdiast;
P_PAm      = 1/3*P_PAsyst + 2/3*P_PAdiast;
P_th       = 0;

% Pulm Veins
P_PUpp    = .2*P_PApp;                            % Pulse pressure in PU (20% of PP_PA) BC REFERENCE
P_PUm     = data.P_PCWave;
P_PUdiast = P_PUm - P_PUpp/3;
P_PUsyst  = P_PUm + 2*P_PUpp/3;

% Systemic Art
P_SAsyst   = data.P_SAsyst;                       % Systolic aortic press (mmHg)    *
P_SAdiast  = data.P_SAdiast;                      % Diastolic aortic press (mmHg)   *
P_SApp     = P_SAsyst - P_SAdiast;
P_SAm      = 1/3*P_SAsyst + 2/3*P_SAdiast;

% Left ventricle
P_LVsyst  = 1.025 * P_SAsyst;                      % Systolic LV pressure (mmHg)
P_LVdiast = 0.975 * P_PUdiast;                     % Diastolic LV pressure (mmHg)

% Syst veins
P_SVdiast  = 1.025*P_RVdiast;
P_SVpp     = .05*P_SApp;                           % Pulse pressure in SV (5% of PP_SA) BC REFERENCE
P_SVsyst   = P_SVdiast + P_SVpp;
P_SVm      = 1/3*P_SVsyst + 2/3*P_SVdiast;

% Other data
Ave_HR     = data.Ave_HR;                         % Average heart rate (beats/min)
gPars.period  = 60/Ave_HR;                        % 1  Period of heart beat (s)
CO         = data.CO_Thermo;                      % Cardiac output Fick (L/min)     *
BW         = data.BW;                             % Body weight (kg)                *
Hgt        = data.Hgt;                            % Height (cm)                     *
Gender     = data.Gender;                         % Gender (1=female, 2=male)       *

if (Gender == 2)
        TotBV = ((0.3669 * (Hgt/100)^3) + (0.03219 * BW) + 0.6041) * 1000;    % Total blood volume (mL)
        BSA   = 0.000579479 * (BW)^(0.38)*(Hgt)^(1.24);                       % Body surface area (m^2)
    else
        TotBV = ((0.3561 * (Hgt/100)^3) + (0.03308 * BW) + 0.1833) * 1000;    % Total blood volume (mL)
        BSA   = 0.000975482 * (BW)^(0.46)*(Hgt)^(1.08);                       % Body surface area (m^2)
end

%  The original Smith model only circulated a portion of the blood so aortic
%  pressure dynamics are not lumped into a general arterila systemic 
%  compartment. Assuming they were simulating a typical 5000 mL total blood 
%  volume they included only 1500 mL (or 30%) in the circulating volume
%  therefore we will multiply our calculated TotBV value by 0.3 to yield 
%  circulating blood volume.

% BC NUMBERS & REFERENCE
V_LVM = 93*BSA - 16;      % Max LV volume (End Diastolic Volume in LV (mL) %BC Reference
SV    = CO/(Ave_HR/60);   % Left stroke volume (mL/beat)
V_LVm = V_LVM - SV;       % Min LV volume (End Systolic Volume in LV mL)

V_RVM = 0.9*V_LVM;        % Max RV volume (10% higher than the  V_LV)  %BC Reference
V_RVm = 0.9*V_LVm;        % Min RV volume (10% higher than the  V_LV)

% Distribution of volume outside the heart 
sart = 0.13;     % Benekin 13
part = 0.03;     % Benekin  3
svein= 0.64;     % Benekin 64  Took volume out of systemic veins and put into the heart
pvein= 0.11;     % Benekin 11 
% Fractions add to 1

% BC REFERENCE
Circ_pa =  0.58 * part*TotBV;    %CircBV; 
Circ_pu =  0.11 * pvein*TotBV;   %CircBV; 
Circ_sa =  0.27 * sart*TotBV;    %CircBV; % was 27 Benekin
Circ_sv =  0.08 * svein*TotBV;   %CircBV; % was 7.5 Benekin

CircBV =  Circ_pa + Circ_pu + Circ_sa + Circ_sv;

% Elastance function driver parameters
H  = Ave_HR;                               % 3  Elastance function parameter (1/s^2)
T2 = gPars.period/2;                       % 4  Elastance function parameter (s)

% Left ventricle free wall parameters
V_d_lvf    = 10;                           % 9  LV ES zero P volume (mL)     c. Smith_Andreassen
E_es_lvf   = P_SAsyst /(V_LVm - V_d_lvf);  % 10 LV free wall elast (mmHg/mL) 
P_0_lvf    = 0.125;                        % 11 LV ED pressure param (mmHg)  c. Smith_Andreassen
lambda_lvf = data.lambda_lvf;              % 12 LV ED pressure param (1/mL)  c. Burkhoff
V_0_lvf    = 10+2;                         % 13 LV ED pressure param (mL)    c. Smith_Andreassen

% Right ventricle free wall parameters
V_d_rvf    = 0.9*V_d_lvf;                  % 14 RV ES zero P volume (mL)     c. Smith_Andreassen
E_es_rvf   = P_PAsyst/(V_RVm -V_d_rvf);    % 15 RV free wall elast (mmHg/mL) 
P_0_rvf    = 0.25;                         % 16 RV ED pressure param (mmHg)  c. Smith_Andreassen
lambda_rvf = data.lambda_rvf;              % 17 RV ED pressure param (1/mL)  c. Burkhoff
V_0_rvf    = 0.9*(10+2);                   % 18 RV ED pressure param (mL)    c. Smith_Andreassen

% Pulmonary artery and vein parameters
E_es_pa = P_PAsyst/Circ_pa;                % 20 Pulm artery elastance (mmHg/mL)  
E_es_pu = P_PUm/Circ_pu;                   % 22 Pulm vein elastance (mmHg/mL)  
R_pul   = (P_PAsyst - P_PUm)/CO;           % 23 Pulm vasc resist (mmHg*s/mL) 

% Systemic artery and vein parameters
E_es_sa = P_SAsyst/Circ_sa;                % 25 Syst art elastance (mmHg/mL)  
E_es_sv = P_SVm/Circ_sv;                   % 27 Syst vein elastance (mmHg/mL)
R_sys   = (P_SAsyst - P_SVm)/CO;           % 28 Syst art resistance (mmHg*s/mL) 

% Heart valve paramenters
R_mt = (P_PUdiast - P_LVdiast)/CO;
R_av = (P_LVsyst  - P_SAsyst)/CO;         
R_tc = (P_SVdiast - P_RVdiast)/CO;
R_pv = (P_RVsyst  - P_PAsyst)/CO;
        
%%%%% INITIAL CONDITIONS %%%%%
V_lv0 = V_LVM - V_d_lvf; 
V_rv0 = V_RVM - V_d_rvf; 
V_pa0 = Circ_pa;
V_pu0 = Circ_pu;
V_sa0 = Circ_sa;
V_sv0 = Circ_sv;
        
pars = [H T2                                      ...  %  1 -  2
E_es_lvf V_d_lvf     P_0_lvf  lambda_lvf  V_0_lvf ...  %  3 -  7
E_es_rvf V_d_rvf     P_0_rvf  lambda_rvf  V_0_rvf ...  %  8 - 12
E_es_pa  E_es_pu  R_pul                           ...  % 13 - 15
E_es_sa  E_es_sv     R_sys                        ...  % 16 - 18
R_mt     R_av        R_tc     R_pv;]';            ...  % 19 - 22

gPars.pars = pars;

pars = log(pars);
hi  = pars + log(4);
low = pars - log(10);

Init = [V_lv0, V_rv0, V_pa0, V_pu0, V_sa0, V_sv0];

end 
