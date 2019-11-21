% ***********************************************************************************
%                       D R I V E R   B A S I C  for
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
% Retrieves nominal values (pars), initial conditions (Init), upper (hi) and lower (low)
% bounds from load_global. 
%
% load res233_opt.mat x: loads the Patient optimized file with optimized 
% parameter values in order to reevaluate the pressures and volumes.
%
% pars = exp(x) >>> exponentiate the optimized parameter set
% 
% FOR COMPUTING NOMINAL PRESSURES AND VOLUME:
% comment out load res233_opt.mat x, uncomment pars = exp(pars). This
% exponentiates the pars from load_global.
% ***********************************************************************************


function DriverBasic()

% Assign the patient data
data = Patient233

% Run load_global to extract initial values (Init) and parameter structure (pars)
[pars,Init,hi,low,gPars] = load_global(data);
pars = exp(pars); %Exponentiate the nominal parameter values

% To load the optimized parameters uncomment the following two lines and
% comment out pars = exp(pars) above.

% load results_opt.mat x
% pars = exp(x);

dt = 0.001;
period       = round(gPars.period/dt)*dt;
gPars.period = period;
gPars.pts    = round(gPars.period/dt);

data.NumBeats =  80;
data.BeatsSS  =  70;
tend = data.NumBeats*gPars.period;
data.time    = [0:dt:tend]; 
data.Periods = [0:gPars.period:data.NumBeats*gPars.period];

data.Init = Init;

% Set ODE/DAE options
options = odeset('AbsTol',gPars.ABS_TOL,'RelTol',gPars.REL_TOL);

% Solve over the steady state time span with ode15s 
sol  = ode15s(@dXdT,[data.time(1) data.time(end)],data.Init,options,pars,gPars);
sols = deval(sol,data.time);

ID  = max(find(data.time    <= data.BeatsSS*gPars.period+1e-6));
IDp = max(find(data.Periods <= data.BeatsSS*gPars.period+1e-6));
timesol = data.time(ID:end);
timeper = data.Periods(IDp:end);

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
P_th        = 0; %- pars(16);                % Thoracic pressure (mmHg)  

% Systemic arteries and vein parameters
E_es_sa     = pars(16);                      % Syst art elastance (kPa/mL)
E_es_sv     = pars(17);                      % Syst vein elastance (kPa/mL)
R_sys       = pars(18);                      % Syst art resistance (kPa*s/mL)
    
% Heart valve paramenters
R_mt        = pars(19);                      % Mitral valve resist (kPa*s/mL)
R_av        = pars(20);                      % Aortic valve resist (kPa*s/mL)
R_tc        = pars(21);                      % Tricuspid vlv resist (kPa*s/mL)
R_pv        = pars(22);                      % Pulmon vlv resist (kPa*s/mL)


%%% P_rv
V_lv        = sols(1,:);
V_rv        = sols(2,:);

tau         = data.time - (floor(data.time/gPars.period) * gPars.period);
e_t         = exp((-1) * H * (tau-T2).^2);

%%% P_rv
P_es_rvf    = E_es_rvf * (V_rv - V_d_rvf);
P_ed_rvf    = P_0_rvf * (exp(lambda_rvf * (V_rv - V_0_rvf)) - 1);
P_rv        = (e_t .* P_es_rvf) + ((1-e_t) .* P_ed_rvf) + P_th;
 
%%% P_lv
P_es_lvf    = E_es_lvf * (V_lv - V_d_lvf);
P_ed_lvf    = P_0_lvf * (exp(lambda_lvf * (V_lv - V_0_lvf)) - 1);
P_lv        = (e_t .* P_es_lvf) + ((1-e_t) .* P_ed_lvf) + P_th;

%%% P_pa
V_pa        = sols(3,:);
P_pa        = E_es_pa * V_pa + P_th;

%%% P_pu
V_pu        = sols(4,:);
P_pu        = E_es_pu * V_pu + P_th;

%%% P_sa
V_sa        = sols(5,:);
P_sa        = E_es_sa * V_sa;

%%% P_sv
V_sv        = sols(6,:);
P_sv        = E_es_sv * V_sv;

%%% CO
Q_sys  = (P_sa - P_sv) / R_sys;
Q_pul  = (P_pa - P_pu) / R_pul;

C1 = data.BeatsSS;
C2 = data.NumBeats;
PPC = gPars.pts;

Periods = data.Periods;
for i = C1:C2
  I1 = (i-1)*PPC+1;
  I2 = i*PPC+1;
  tdc   = data.time(I1:I2);
  P_sac = P_sa(I1:I2);
  P_pac = P_pa(I1:I2);
  P_rvc = P_rv(I1:I2);
  V_lvc = V_lv(I1:I2);
  V_rvc = V_rv(I1:I2);
  
  PsaSys(i-C1+1) = max(P_sac); % systemic arteries
  PsaDia(i-C1+1) = min(P_sac);
  PpaSys(i-C1+1) = max(P_pac); % pulmonary arteries
  PpaDia(i-C1+1) = min(P_pac);
  PrvSys(i-C1+1) = max(P_rvc); % Right ventricle
  PrvDia(i-C1+1) = min(P_rvc);
  CO_sys(i-C1+1) = trapz(tdc,Q_sys(I1:I2))/(tdc(end)-tdc(1));
  CO_pul(i-C1+1) = trapz(tdc,Q_pul(I1:I2))/(tdc(end)-tdc(1));
  
  CO_Vl = (max(V_lvc) - min(V_lvc))/(tdc(end)-tdc(1));
  CO_Vr = (max(V_rvc) - min(V_rvc))/(tdc(end)-tdc(1));
  
  display([tdc(1) tdc(end) CO_sys(i-C1+1) CO_pul(i-C1+1) CO_Vl CO_Vr])
  Ppum(i-C1+1)   = trapz(tdc,P_pu(I1:I2))/(tdc(end)-tdc(1));  
end;

% Data
P_SAsys   = data.P_SAsyst *ones(size(PsaSys));   % Systolic aortic press (mmHg)    
P_SAdia   = data.P_SAdiast*ones(size(PsaDia));   % Diastolic aortic press (mmHg)   
P_PAsys   = data.P_PAsyst *ones(size(PpaSys));   % Systolic PA press (mmHg)    
P_PAdia   = data.P_PAdiast*ones(size(PpaDia));   % Diastolic PA press (mmHg)   
P_RVsys   = data.P_RVsyst *ones(size(PrvSys));   % Systolic Right Ventricle press (mmHg)   
P_RVdia   = data.P_RVdiast*ones(size(PrvDia));   % Diastolic Right Ventricle  press (mmHg)  
COd       = data.CO_Thermo*ones(size(CO_sys));
PCWp      = data.P_PCWave*ones(size(Ppum));


Vtotal = V_lv + V_rv + V_pa + V_pu + V_sa + V_sv;
fn = strcat(results,'.mat'); %If you are running with optimized results, 
% change 'results' to  'results_opt'
save(fn);

% Figure(1) plots 
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
