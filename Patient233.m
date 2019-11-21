function data = Patient233
% Patient 233 with right heart catheterization on 29 October 2015
P_RVsyst  = 40;                             % Systolic RV pressure (mmHg)     *
P_RVdiast = 4;                              % Diastolic RV pressure (mmHg)    *
P_PAsyst  = 35;                             % Systolic pulm art press (mmHg)  *
P_PAdiast = 19;                             % Diast pulm art press (mmHg)     *
P_PCWave  = 0.98*19;                        % Average pulm wedge press (mmHg) * 
P_SAsyst  = 149;                            % Systolic arterial press (mmHg)    *
P_SAdiast = 83;                             % Diastolic arterial press (mmHg)   *
Ave_HR    = 83.33;                          % Average heart rate (beats/min)
CO_Fick   = 7.34*1000/60;                   % Cardiac output Fick (mL/sec)     *
CO_Thermo = 10.33*1000/60;                  % Crdiac outpt thermodil (mL/sec)  *
BW        = 96.1;                           % Body weight (kg)                *
Hgt       = 172;                            % Height (cm)                     *
Gender    = 1;                              % Gender (1=female, 2=male)       *

lambda_lvf = 0.03;                          % 0.03 c. Burkhoff but some patients are adjusted to fit the data
lambda_rvf = 0.025;                         % 0.025 c. Burkhoff but some patients are adjusted to fit the data

Patient233_Values = {P_RVsyst P_RVdiast P_PAsyst P_PAdiast P_PCWave P_SAsyst P_SAdiast ...
    Ave_HR CO_Fick CO_Thermo BW Hgt Gender lambda_lvf lambda_rvf};

Patient233_Fields = {'P_RVsyst' 'P_RVdiast' 'P_PAsyst' 'P_PAdiast' 'P_PCWave' 'P_SAsyst' 'P_SAdiast' ...
    'Ave_HR' 'CO_Fick' 'CO_Thermo' 'BW' 'Hgt' 'Gender' 'lambda_lvf' 'lambda_rvf'};
  
data = cell2struct(Patient233_Values, Patient233_Fields, 2);
clear Patient233_Values Patient233_Fields
end
