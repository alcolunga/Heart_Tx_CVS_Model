% ***********************************************************************************
%                    d X d T   O D E   F U N C T I O N   for
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
%   This function contains the algebraic and differential expressions that describe
%
%   Model originally created on     17 January 2016
%   Model last modfied on           24 October 2016
%
%   Reproduced by       Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
% ***********************************************************************************

function dV = dXdT(time,y,pars,gPars)
   
    % Elastance function driver parameters
    H           = pars(1);                       % Elastance fctn param (1/s^2)
    T2          = pars(2);                       % Elastance fctn param (s)
    
    % Left ventricle free wall parameters
    E_es_lvf    = pars(3);                       % LV free wall elastance (kPa/mL)
    V_d_lvf     = pars(4);                       % LV ES zero P volume (mL)
    P_0_lvf     = pars(5);                       % LV ED pressure param (kPa)
    lambda_lvf  = pars(6);                       % LV ED pressure param (1/mL)
    V_0_lvf     = pars(7);                       % LV ED pressure param (mL)
    
    % Right ventricle free wall parameters
    E_es_rvf    = pars(8);                       % RV free wall elastance (kPa/mL)
    V_d_rvf     = pars(9);                       % RV ES zero P volume (mL)
    P_0_rvf     = pars(10);                      % RV ED pressure param (kPa)
    lambda_rvf  = pars(11);                      % RV ED pressure param (1/mL)
    V_0_rvf     = pars(12);                      % RV ED pressure param (mL)
    
    % Pulmonary artery and vein parameters
    E_es_pa     = pars(13);                      % Pulm artery elastance (kPa/mL)
    E_es_pu     = pars(14);                      % Pulm vein elastance (kPa/mL)
    R_pul       = pars(15);                      % Pulm vasc resistance (kPa*s/mL)
    P_th        = 0; %-pars(16);
    
    % Aortic and vena cava parameters
    E_es_sa     = pars(16);                      % Syst Art elastance (kPa/mL)
    E_es_sv     = pars(17);                      % Syst Vein elastance (kPa/mL)
    R_sys       = pars(18);                      % Syst art resistance (kPa*s/mL)
    
    % Heart valve paramenters
    R_mt        = pars(19);                      % Mitral valve resist (kPa*s/mL)
    R_av        = pars(20);                      % Aortic valve resist (kPa*s/mL)
    R_tc        = pars(21);                      % Tricuspid vlv resist (kPa*s/mL)
    R_pv        = pars(22);                      % Pulmon vlv resist (kPa*s/mL)
    
     
    % Unpack the X vector
    V_lv = y(1);
    V_rv = y(2);
    V_pa = y(3);
    V_pu = y(4);
    V_sa = y(5);
    V_sv = y(6);
    
    %   <component name="driver_function">
    tau = time - (floor(time/gPars.period) * gPars.period);
    e_t = exp((-1) * H * (tau-T2)^2);

    %   <component name="lvf_calculator">
    P_es_lvf = E_es_lvf * (V_lv - V_d_lvf);
    P_ed_lvf = P_0_lvf * (exp(lambda_lvf * (V_lv - V_0_lvf)) - 1);
    P_lv     = (e_t * P_es_lvf) + ((1-e_t) * P_ed_lvf) + P_th;
    
    %   <component name="rvf_calculator">
    P_es_rvf = E_es_rvf * (V_rv - V_d_rvf);
    P_ed_rvf = P_0_rvf * (exp(lambda_rvf * (V_rv - V_0_rvf)) - 1);
    P_rv     = (e_t * P_es_rvf) + ((1-e_t) * P_ed_rvf) + P_th;
    
    %   <component name="pulmonary_artery and vein">
    P_pa = E_es_pa * (V_pa) + P_th;
    P_pu = E_es_pu * (V_pu) + P_th;
    Q_pul = (P_pa-P_pu) / R_pul;

    %   <component name="syst art and syst vein">
    P_sa = E_es_sa * (V_sa);
    P_sv = E_es_sv * (V_sv);
    Q_sys = (P_sa-P_sv) / R_sys;
    if P_pu > P_lv
        Q_mt = (P_pu - P_lv)/R_mt;
    else
        Q_mt = 0;
    end
    if P_lv > P_sa
        Q_av = (P_lv - P_sa)/R_av;
    else
        Q_av = 0;
    end
    if P_sv > P_rv
        Q_tc = (P_sv - P_rv)/R_tc;
    else
        Q_tc = 0;
    end
    if P_rv > P_pa
        Q_pv = (P_rv - P_pa)/R_pv;
    else
        Q_pv = 0;
    end
    dVlv = Q_mt  - Q_av;
    dVrv = Q_tc  - Q_pv;
    dVpa = Q_pv  - Q_pul;
    dVpu = Q_pul - Q_mt;
    dVsa = Q_av  - Q_sys;
    dVsv = Q_sys - Q_tc;

    dV = [dVlv dVrv dVpa dVpu dVsa dVsv]';     
end

