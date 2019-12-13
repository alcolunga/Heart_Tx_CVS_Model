% ***********************************************************************************
%                       D R I V E R   B A S I C  for
%            C A R D I O V A S C U L A R   S Y S T E M S   M O D E L
% ***********************************************************************************
% Copyright (C) 2019 A.L. Colunga(1), N.P. Woodall(1), M.S. Olufsen(1), and
% B.E. Carlson(2)
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
% 
% ***********************************************************************************
% load res233_opt.mat xopt: loads the Patient optimized file with optimized 
% parameter values in order to reevaluate the pressures and volumes.
%
% INDMAP  = [6 11 13 14 15 16 18 22] the parameters optimized
%           (6 lambda_lvf, 11 lambda_rvf, 13 E_pa, 14 E_pu, 15 R_pul,
%            16 E_sa, 18 R_sys, 22 R_pv)
%
% pars(INDMAP) = exp(xopt) >>> exponentiate the parameter set
%
% FOR COMPUTING NOMINAL PRESSURES AND VOLUME:
% comment out load results_opt.mat xopt
%
% Input:
% none
% 
% Output: 
% Results       mat file 
%
% Dependencies: 
% load_global   write inputs and outputs in sentence 
% dXdT          same here as above
% ***********************************************************************************

function DriverBasic()

% Assign the patient data
data = Patient233;

% Compute nominal parameter values and inital conditions
[gPars,Init] = load_global(data);
pars = exp(gPars.pars); %Exponentiate the log-scaled nominal parameter values

% %Optimized parameters:(6) lambda_lvf, (11) lambda_rvf, (13) E_pa, (14) E_pu, 
% %(15) R_pul, (16) E_sa, (18) R_sys, (22) R_pv.
%
% % Load optimized parameters and assign their index value for the entire
% $ parameter set and exponentiate.
% % To load uncomment the following three lines

%load results_opt.mat xopt 
%INDMAP  = [6 11 13 14 15 16 18 22];
%pars(INDMAP) = exp(xopt)

% Set up simulation time 
dt = 0.001;                                % Solution is displayed at time-step dt = 0.001
period       = round(gPars.period/dt)*dt;  % Length of cardiac cycle

% 
gPars.period = period;                     % Cardiac cycle lenght added to data vector
gPars.pts    = round(gPars.period/dt);     % Number of time-points per cardiac cycle

% Generating a structure 
data.NumBeats =  80;                       % Number of cycles computed
data.BeatsSS  =  70;                       % Results displayed after steady state has been obtained
tend = data.NumBeats*gPars.period;         % Last time point
data.time    = [0:dt:tend];                % Time vector
data.Periods = [0:gPars.period:data.NumBeats*gPars.period]; % Vector with times for periods

data.Init = Init;                          % Adding initial conditions to data structure

% Options for ODE solver
options = odeset('AbsTol',gPars.ABS_TOL,'RelTol',gPars.REL_TOL);

% Solve ODEs over the steady state time span with ode15s 
sol  = ode15s(@dXdT,[data.time(1) data.time(end)],data.Init,options,pars,gPars);
sols = deval(sol,data.time);

% Generate time vecotr for displaying solution starting at t=data.BeatsSS (70s)
ID  = max(find(data.time    <= data.BeatsSS*gPars.period+1e-6));
IDp = max(find(data.Periods <= data.BeatsSS*gPars.period+1e-6));
timesol = data.time(ID:end);
timeper = data.Periods(IDp:end);

% Elastance function driver parameters
H           = pars(1);                       % Elastance function (1/s^2)
T2          = pars(2);                       % Elastance function (s)
    
% Left ventricle free wall parameters
E_lvf       = pars(3);                       % LV free wall elastance (kPa/mL)
V_d_lvf     = pars(4);                       % LV end systolic zero pressure volume (mL)
P_0_lvf     = pars(5);                       % LV end diastolic pressure (kPa)
lambda_lvf  = pars(6);                       % LV end diastolic pressure (1/mL)
V_0_lvf     = pars(7);                       % LV ebd diastolic pressure (mL)
    
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

%Thorax
P_th        = 0;                             % Thoracic pressure (mmHg)  

% Systemic arteries and vein parameters
E_sa        = pars(16);                      % SA elastance (kPa/mL)
E_sv        = pars(17);                      % SV elastance (kPa/mL)
R_sys       = pars(18);                      % Systolic vascular resistance (kPa*s/mL)
    
% Heart valve paramenters
R_mt        = pars(19);                      % Mitral valve resistance (kPa*s/mL)
R_av        = pars(20);                      % Aortic valve resistance (kPa*s/mL)
R_tc        = pars(21);                      % Tricuspid valve resistance (kPa*s/mL)
R_pv        = pars(22);                      % Pulmon valve resistance (kPa*s/mL)


%%% Ventricles 
V_lv        = sols(1,:);                     % Left ventricular volume (mL)
V_rv        = sols(2,:);                     % Right ventricular volume (mL)

tau         = data.time - (floor(data.time/gPars.period) * gPars.period);
e_t         = exp((-1) * H * (tau-T2).^2);   % Elastance function

%%% Right ventricular pressure
P_es_rvf       = E_rvf * (V_rv - V_d_rvf);
P_ed_rvf       = P_0_rvf * (exp(lambda_rvf * (V_rv - V_0_rvf)) - 1);
P_rv        = (e_t .* P_es_rvf) + ((1-e_t) .* P_ed_rvf) + P_th;
 
%%% Left ventricular pressure
P_es_lvf       = E_lvf * (V_lv - V_d_lvf);
P_ed_lvf       = P_0_lvf * (exp(lambda_lvf * (V_lv - V_0_lvf)) - 1);
P_lv        = (e_t .* P_es_lvf) + ((1-e_t) .* P_ed_lvf) + P_th;

%%% Pulmonary arteries (PA)
V_pa        = sols(3,:);                    % Volume (mL)
P_pa        = E_pa * V_pa + P_th;           % Pressure (mmHg)

%%% Pulmonary veins (PU)
V_pu        = sols(4,:);                    % Volume (mL)
P_pu        = E_pu * V_pu + P_th;           % Pressure (mmHg)

%%% Systemic arteries (SA)
V_sa        = sols(5,:);                    % Volume (mL)
P_sa        = E_sa * V_sa;                  % Pressure (mmHg)

%%% Systemic veins (SV)
V_sv        = sols(6,:);                    % Volume (mL)
P_sv        = E_sv * V_sv;                  % Pressure (mmHg)

%%% Cardiac output (CO)
Q_sys  = (P_sa - P_sv) / R_sys;             % Systemic cardiac output (mL/s)
Q_pul  = (P_pa - P_pu) / R_pul;             % Pulmonary cardiac output (mL/s)

%%% Information extracted over each cardiac cycle displayed (from t=70 to 80)
C1 = data.BeatsSS;
C2 = data.NumBeats;
PPC = gPars.pts;

Periods = data.Periods;
for i = C1:C2  
  I1 = (i-1)*PPC+1;
  I2 = i*PPC+1;
  tdc   = data.time(I1:I2);                 % Time points in current cycle
  P_sac = P_sa(I1:I2);                      % SA pressure in current cycle
  P_pac = P_pa(I1:I2);                      % PA pressure in current cycle
  P_rvc = P_rv(I1:I2);                      % RV pressure in current cycle
  V_lvc = V_lv(I1:I2);                      % LV volume in current cycle
  V_rvc = V_rv(I1:I2);                      % RV volume in current cycle
  
  PsaSys(i-C1+1) = max(P_sac);              % SA systolic pressure
  PsaDia(i-C1+1) = min(P_sac);              % SA diastolic pressure
  PpaSys(i-C1+1) = max(P_pac);              % PA systolic pressure
  PpaDia(i-C1+1) = min(P_pac);              % PA diastolic pressure
  PrvSys(i-C1+1) = max(P_rvc);              % Right ventriclular systolic pressure
  PrvDia(i-C1+1) = min(P_rvc);              % Right ventriclular diastolic pressure
  CO_sys(i-C1+1) = trapz(tdc,Q_sys(I1:I2))/(tdc(end)-tdc(1)); % Systemic cardiac output
  CO_pul(i-C1+1) = trapz(tdc,Q_pul(I1:I2))/(tdc(end)-tdc(1)); % Pulmonary cardiac output
  
  CO_Vl = (max(V_lvc) - min(V_lvc))/(tdc(end)-tdc(1)); % Cardiac output left ventricle
  CO_Vr = (max(V_rvc) - min(V_rvc))/(tdc(end)-tdc(1)); % Cardiac output right ventricle
  % Note all calculations of cardiac output should be the same.
  
  display([tdc(1) tdc(end) CO_sys(i-C1+1) CO_pul(i-C1+1) CO_Vl CO_Vr]) % Display values for current cardiac cycle
  Ppum(i-C1+1) = trapz(tdc,P_pu(I1:I2))/(tdc(end)-tdc(1));  % ??????
end;

% Data
P_SAsys   = data.P_SAsyst *ones(size(PsaSys));   % Systolic aortic press (mmHg)    
P_SAdia   = data.P_SAdiast*ones(size(PsaDia));   % Diastolic aortic press (mmHg)   
P_PAsys   = data.P_PAsyst *ones(size(PpaSys));   % Systolic PA press (mmHg)    
P_PAdia   = data.P_PAdiast*ones(size(PpaDia));   % Diastolic PA press (mmHg)   
P_RVsys   = data.P_RVsyst *ones(size(PrvSys));   % Systolic Right Ventricle press (mmHg)   
P_RVdia   = data.P_RVdiast*ones(size(PrvDia));   % Diastolic Right Ventricle  press (mmHg)  
COd       = data.CO_Thermo*ones(size(CO_sys));   % Cardiac output (mL/s)
PCWp      = data.P_PCWave*ones(size(Ppum));      % Pulmonary capilary wedge pressure (mmHg)


% Total stressed volume (should be constant)
Vtotal = V_lv + V_rv + V_pa + V_pu + V_sa + V_sv;

% Save results to file
save Results233.mat;

% Figure(1) 
% Left ventricle, pulmonary vein, & aorta comparison of computed 
% results (solid lines) and data (broken lines) (2,3,1) 
%
% Right ventricle, pulmonary artery, & vena cava comparison of 
% computed results (solid lines) and data (broken lines) (2,3,2) 
%
% comparison of computed cardiac output and data. (2,3,3)
%
% comparison of left and right ventricular pressure. (2,3,4)
%
% computed left and right ventricular volume (no data available).(2,3,5)
%
% Left and right ventricular pressure volume loop (PV loop) with calculated 
% stroke work. (2,3,6)

t1 = timesol(1);
t2 = timesol(1)+5;
figure(1);
subplot(2,3,1);hold on;
    h=plot(timesol,P_lv(ID:end),timesol,P_sa(ID:end),timesol,P_pu(ID:end)); hold on;
    set(h(1),'Color',[.7 .06 .06])  
    set(h(2),'Color',[1 .50 .50])
    set(h(3),'Color',[1 .54 0])
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Time (sec)');
    ylabel('Pressure (mmHg)');
    g=plot(timeper,PCWp,'--','linewidth',3);
    set(g,'Color',[1 .54 0])
    g=plot(timeper,P_SAsys,'--',timeper,P_SAdia,'--','linewidth',3);
    set(g,'Color',[1 .50 .50])
    legend('P_{lv}','P_{sa}','P_{pv}');
    grid on;
    xlim([t1 t2]);    
    ylim([0 175]);
    
subplot(2,3,2);hold on;
    h = plot(timesol,P_rv(ID:end),'b',timesol,P_pa(ID:end),'c',timesol,P_sv(ID:end),'g'); hold on
    g = plot(timeper,P_PAsys,'c--',timeper,P_PAdia,'c--');hold on;
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Time (sec)');
    ylabel('Pressure (mmHg)');
    plot(timeper,P_RVsys,'b--',timeper,P_RVdia,'b--','linewidth',3);
    legend('P_{rv}','P_{pa}','P_{sv}');
    xlim([t1 t2]);  
    ylim([0 50]);
    grid on;
    
subplot(2,3,3);hold on;
    g=plot(timeper,COd,'--k',timeper,CO_sys,'--r',timeper,CO_pul,'--b','linewidth',3);
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Time (sec)');
    ylabel('Flow (mL/sec)');
    legend('CO_{data}','CO_{sys}','CO_{pul}');
    grid on;
    xlim([t1 t2]);
   
    
subplot(2,3,4);hold on;
    h=plot(timesol,P_lv(ID:end),timesol,P_rv(ID:end),'b');
    set(h(1),'Color',[.7 .06 .06])
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Time (sec)');
    ylabel('Pressure (mmHg)');
    legend('P_{lv}','P_{rv}');
    grid on;
    xlim([t1 t2]);
    ylim([0 165]);
    
subplot(2,3,5);hold on;
    h=plot(timesol,V_lv(ID:end),timesol,V_rv(ID:end),'b');
    set(h(1),'Color',[.7 .06 .06])
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Time (sec)');
    ylabel('Volume (mL)');
    legend('V_{lv}','V_{rv}');
    grid on;
    xlim([t1 t2]);
    ylim([40 165]);

subplot(2,3,6);hold on;
    h=plot(V_rv(ID:end),P_rv(ID:end),'b');
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    grid on;
    hold on;
    h=plot(V_lv(ID:end),P_lv(ID:end));
    set(h,'Color',[.7 .06 .06])
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Volume (mL)');
    ylabel('Pressure (mmHg)');
    grid on;
    legend('RV','LV')
    ylim([0 167]);
    
end
