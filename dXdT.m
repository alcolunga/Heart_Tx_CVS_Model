% ***********************************************************************************
%                    d X d T   O D E   F U N C T I O N   for
%            C A R D I O V A S C U L A R   S Y S T E M S   M O D E L
% ***********************************************************************************
% Copyright (C) 2019 A.L. Colunga(1), N.P. Woodall(1), M.S. Olufsen(1), and B.E. Carlson
% 
% 1) Department of Mathematics, NC State University 
% 2) Molecular and Integrative Physiology, University of Michigan
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software"), to deal in
% the Software without restriction, including without limitation the rights to
% use, copy, modify, and merge the Software subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% Any work resulting from the use of the Software should cite the associated 
% article arXiv:1812.11857
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, OR DATA, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
% ***********************************************************************************

% ***********************************************************************************
%   This function contains the algebraic and differential expressions that describe
% CHANGE BELOW
% Input:
% data   structure with patient data generated by Patient233.m
% 
% Output: 
% gPars  structure with nominal parameter values
% Init   vector with initial conditions
%
% References:
% Benekin 
% Smith & Andreassen
% Burkhoff
% ***********************************************************************************

function dV = dXdT(time,y,pars,gPars)
   
    % Elastance function driver parameters
    H           = pars(1);                       % Elastance fctn param (1/s^2)
    T2          = pars(2);                       % Elastance fctn param (s)
    
    % Left ventricle free wall parameters
    E_lvf       = pars(3);                       % LV free wall elastance (kPa/mL)
    V_d_lvf     = pars(4);                       % LV end systolic zero pressure volume (mL)
    P_0_lvf     = pars(5);                       % LV end diastolic pressure (kPa)
    lambda_lvf  = pars(6);                       % LV end diastolic pressure (1/mL)
    V_0_lvf     = pars(7);                       % LV end diastolic pressure (mL)
    
    % Right ventricle free wall parameters
    E_rvf       = pars(8);                       % RV free wall elastance (kPa/mL)
    V_d_rvf     = pars(9);                       % RV end systolic zero pressure volume (mL)
    P_0_rvf     = pars(10);                      % RV end diastolic pressure (kPa)
    lambda_rvf  = pars(11);                      % RV end diastolic pressure (1/mL)
    V_0_rvf     = pars(12);                      % RV end diastolic pressure (mL)
    
    % Pulmonary artery and vein parameters
    E_pa        = pars(13);                      % PA elastance (kPa/mL)
    E_pu        = pars(14);                      % PU elastance (kPa/mL)
    R_pul       = pars(15);                      % Pulmonary vascular resistance (kPa*s/mL)
    P_th        = 0; %-pars(16);
    
    % Aortic and vena cava parameters
    E_sa        = pars(16);                      % SA elastance (kPa/mL)
    E_sv        = pars(17);                      % SV elastance (kPa/mL)
    R_sys       = pars(18);                      % Systemic vascular resistance (kPa*s/mL)
    
    % Heart valve paramenters
    R_mt        = pars(19);                      % Mitral valve resist (kPa*s/mL)
    R_av        = pars(20);                      % Aortic valve resist (kPa*s/mL)
    R_tc        = pars(21);                      % Tricuspid vlv resist (kPa*s/mL)
    R_pv        = pars(22);                      % Pulmon vlv resist (kPa*s/mL)
    
     
    % Unpack the state variable vector
    V_lv = y(1);
    V_rv = y(2);
    V_pa = y(3);
    V_pu = y(4);
    V_sa = y(5);
    V_sv = y(6);
    
    % Elastance function 
    tau = time - (floor(time/gPars.period) * gPars.period);         % Timing for ventricular contraction
    e_t = exp((-1) * H * (tau-T2)^2);                               % Elastance function

    % Left ventricular pressure
    P_es_lvf  = E_lvf * (V_lv - V_d_lvf);                            % End systolic pressure (mmHg)
    P_ed_lvf  = P_0_lvf * (exp(lambda_lvf * (V_lv - V_0_lvf)) - 1);  % End diastolic pressure (mmHg)
    P_lv   = (e_t * P_es_lvf) + ((1-e_t) * P_ed_lvf) + P_th;            % LV pressure (mmHg)
    
    % Right ventricular pressure
    P_es_rvf = E_rvf * (V_rv - V_d_rvf);                         % End systolic pressure (mmHg)
    P_ed_rvf = P_0_rvf * (exp(lambda_rvf * (V_rv - V_0_rvf)) - 1);  % End diastolic pressure (mmHg)
    P_rv     = (e_t * P_es_rvf) + ((1-e_t) * P_ed_rvf) + P_th;      % RV pressure (mmHg)
    
    % Pulmonary pressure and flow
    P_pa = E_pa * (V_pa) + P_th;     % PA pressure (mmHg)
    P_pu = E_pu * (V_pu) + P_th;     % PU pressure (mmHg)
    Q_pul = (P_pa-P_pu) / R_pul;        % Pulmonary flow (mL/s)

    % Systemic pressure and flow
    P_sa = E_sa * (V_sa);            % SA pressure (mmHg)
    P_sv = E_sv * (V_sv);            % SV pressure (mmHg)
    Q_sys = (P_sa-P_sv) / R_sys;        % Systemic flow (mL/s)
    
    % Valve flow
    if P_pu > P_lv                      % Flow through mitral valve (mL/s)
        Q_mt = (P_pu - P_lv)/R_mt;
    else
        Q_mt = 0;
    end
    if P_lv > P_sa                      % Flow through aortic valve (mL/s)
        Q_av = (P_lv - P_sa)/R_av;
    else
        Q_av = 0;
    end
    if P_sv > P_rv                      % Flow through tricuspid valve (mL/s)
        Q_tc = (P_sv - P_rv)/R_tc;
    else
        Q_tc = 0;
    end
    if P_rv > P_pa                      % Flow through pulmonary valve (mL/s)
        Q_pv = (P_rv - P_pa)/R_pv;
    else
        Q_pv = 0;
    end
    
    % Change in volumes
    dVlv = Q_mt  - Q_av;           % dLV volume 
    dVrv = Q_tc  - Q_pv;           % dRV volume
    dVpa = Q_pv  - Q_pul;          % dPA volume
    dVpu = Q_pul - Q_mt;           % dPU volume
    dVsa = Q_av  - Q_sys;          % dSA volume
    dVsv = Q_sys - Q_tc;           % dSV volume

    % Vector returned to ODE solver
    dV = [dVlv dVrv dVpa dVpu dVsa dVsv]';     
end

