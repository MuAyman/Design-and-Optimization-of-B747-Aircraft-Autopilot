clear; clc; close all;

%% Reading Aircraft Data From Excel
Aircraft_data = xlsread('B-747_data.xlsx','B1:B50');

%% Flight Conditoins Data
flight = num2cell(Aircraft_data(1:18));
[W,g,s,b,c,h,M,V_tot,CG,Ix,Iy,Iz,Ixz,Ixy,Iyz,Alpha,q,u]= deal(flight{:});
m = W/g;
%% Lateral Dinemntional Stability Derivatives Data
lateral = num2cell(Aircraft_data(36:50));
[Yv,Lv,Nv,Yp,Lp,Np,Yr,Lr,Nr,YCA,LCA,NCA,YCR,LCR,NCR] = deal(lateral{:});

%% Lateral non-Dinemntional Stability Derivatives Data
Cyv = Yv*m*u/(q*s);
Clv = Lv*Ix*u/(q*s*b);
Cnv = Nv*Ix*u/(q*s*b);
Cyp = 2*Yp*m*u/(q*s*b);
Clp = 2*Lp*Ix*u/(q*s*b^2);
Cnp = 2*Np*Ix*u/(q*s*b^2);
Cyr = 2*Yr*m*u/(q*s*b);
Clr = 2*Lr*Ix*u/(q*s*b^2);
Cnr = 2*Nr*Ix*u/(q*s*b^2);
Cyca = YCA*m/(q*s);
Clca = LCA*Ix/(q*s*b);
Cnca = NCA*Ix/(q*s*b);
Cycr = YCR*m/(q*s);
Clcr = LCR*Ix/(q*s*b);
Cncr = NCR*Ix/(q*s*b);

%% Lateral State Space Model
ALateralFull = [Yv Yp Yr-u g 0; Lv Lp Lr 0 0; Nv Np Nr 0 0; 0 1 0 0 0; 0 0 1 0 0];
BLateralFull = [YCA YCR;LCA LCR;NCA NCR; 0 0; 0 0];
CLateralFull = eye(5);
DLateralFull = zeros(5,2);
LateralSSFull = ss(ALateralFull,BLateralFull,CLateralFull,DLateralFull);
LateralSSFull.OutputName = {'v';'p';'r';'phi';'psi'};
LateralSSFull.InputName = {'Da';'Dr'};

%% Lateral Dutch Roll, Rolling & Spiral modes
% Dutch Roll Mode
ALateralDutch = [Yv Yp Yr-u; Lv Lp Lr;Nv Np Nr];
BLateralDutch = [YCA YCR;LCA LCR;NCA NCR];
CLateralDutch = eye(3);
DLateralDutch = zeros(3,2);
LateralSSDutch = ss(ALateralDutch,BLateralDutch,CLateralDutch,DLateralDutch);
LateralSSDutch.OutputName = {'v';'p';'r'};
LateralSSDutch.InputName = {'Da';'Dr'};
% Rolling Mode
ALateralRolling = Lp;
BLateralRolling = LCA;
CLateralRolling = 1;
DLateralRolling = 0;
LateralSSRolling = ss(ALateralRolling,BLateralRolling,CLateralRolling,DLateralRolling);
LateralSSRolling.OutputName = {'p'};
LateralSSRolling.InputName = {'Da'};
% Spiral Mode
ALateralSpiral = [Lp Lr 0;Np Nr 0; 1 0 0];
BLateralSpiral = [LCR; NCR;0];
CLateralSpiral = eye(3);
DLateralSpiral = zeros(3,1);
LateralSSSpiral = ss(ALateralSpiral,BLateralSpiral,CLateralSpiral,DLateralSpiral);
LateralSSSpiral.OutputName = {'p';'r';'phi'};
LateralSSSpiral.InputName = {'Dr'};

%% Pole Placement In Case Of Insatiability
PolesFull = pole(LateralSSFull);
if ~(isreal(sqrt(-real(PolesFull))))  % enters iff there is a positive pole
    % Full
    PolesFullNew = complex(-abs(real(PolesFull)),imag(PolesFull));
    kFull = place(ALateralFull,BLateralFull,PolesFullNew);
    ALateralFull = ALateralFull-BLateralFull*kFull;
    LateralSSFull = ss(ALateralFull,BLateralFull,CLateralFull,DLateralFull);
    PolesFull = pole(LateralSSFull);
    LateralSSFull.OutputName = {'v';'p';'r';'phi';'psi'};
    LateralSSFull.InputName = {'Da';'Dr'};
    % Dutch Roll Mode
    PolesDutch = pole(LateralSSDutch);
    PolesDutchNew = complex(-abs(real(PolesDutch)),imag(PolesDutch));
    kDutch = place(ALateralDutch,BLateralDutch,PolesDutchNew);
    ALateralDutch = ALateralDutch-BLateralDutch*kDutch;
    LateralSSDutch = ss(ALateralDutch,BLateralDutch,CLateralDutch,DLateralDutch);
    LateralSSDutch.OutputName = {'v';'p';'r'};
    LateralSSDutch.InputName = {'Da';'Dr'};
    % Rolling Mode
    PolesRolling = pole(LateralSSRolling);
    PolesRollingNew = complex(-abs(real(PolesRolling)),imag(PolesRolling));
    kRolling = place(ALateralRolling,BLateralRolling,PolesRollingNew);
    ALateralRolling = ALateralRolling-BLateralRolling*kRolling;
    LateralSSRolling = ss(ALateralRolling,BLateralRolling,CLateralRolling,DLateralRolling);
    LateralSSRolling.OutputName = {'p'};
    LateralSSRolling.InputName = {'Da'};
    % Spiral Mode
    PolesSpiral = pole(LateralSSSpiral);
    PolesSpiralNew = complex(-abs(real(PolesSpiral)),imag(PolesSpiral));
    kSpiral = place(ALateralSpiral,BLateralSpiral,PolesSpiralNew);
    ALateralSpiral = ALateralSpiral-BLateralSpiral*kSpiral;
    LateralSSSpiral = ss(ALateralSpiral,BLateralSpiral,CLateralSpiral,DLateralSpiral);
    LateralSSSpiral.OutputName = {'p';'r';'phi'};
    LateralSSSpiral.InputName = {'Dr'};
end

%% Lateral Eigenvalues, eigenvectors
[eigVectorsLateralFull, eigValuesLateralFull] = eig(ALateralFull);

%% Full & Modes Characteristics (undamped natural frequency, damping ratio, settling time, etc.)
% Full
[WnFull,ZetaFull] = damp(LateralSSFull);
StepInfoFull = stepinfo(step(LateralSSFull));
% Dutch Roll Mode
[WnDutch,ZetaDutch] = damp(LateralSSDutch);
StepInfoDutch = stepinfo(step(LateralSSDutch));
% Rolling Mode
[WnRolling,ZetaRolling] = damp(LateralSSRolling);
StepInfoRolling = stepinfo(step(LateralSSRolling));
% Spiral Mode
[WnSpiral,ZetaSpiral] = damp(LateralSSSpiral);
StepInfoSpiral = stepinfo(step(LateralSSSpiral));

%% Modes Approx. Error
ErorrDutch = abs(WnDutch-WnFull(3:5))/WnFull(3:5)*100;
ErorrRolling = abs(WnRolling-WnFull(5))/WnFull(5)*100;
ErorrSpiral = abs(WnSpiral-WnFull([1 2 5]))/WnFull([1 2 5])*100;

%% Lateral Transfer Function
LateraltfFull = tf(LateralSSFull);
LateraltfFull_v_da = LateraltfFull(1,1);
LateraltfFull_v_dr = LateraltfFull(1,2);
LateraltfFull_p_da = LateraltfFull(2,1);
LateraltfFull_p_dr = LateraltfFull(2,2);
LateraltfFull_r_da = LateraltfFull(3,1);
LateraltfFull_r_dr = LateraltfFull(3,2);
LateraltfFull_phi_da = LateraltfFull(4,1);
LateraltfFull_phi_dr = LateraltfFull(4,2);
LateraltfFull_psi_da = LateraltfFull(5,1);
LateraltfFull_psi_dr = LateraltfFull(5,2);

%% Lateral Modes Transfer Function
% Dutch Roll Mode
LateraltfDutch = tf(LateralSSDutch);
LateraltfDutch_v_da = LateraltfDutch(1,1);
LateraltfDutch_v_dr = LateraltfDutch(1,2);
LateraltfDutch_r_da = LateraltfDutch(2,1);
LateraltfDutch_r_dr = LateraltfDutch(2,2);
% Rolling Mode
LongtfRolling = tf(LateralSSRolling);
% Spiral Mode
LongtfSpiral = tf(LateralSSSpiral);


%% Simulating the system
% Full Step
step(LateralSSFull);
title('Full sys Step Response')
grid on
% Full Impulse
figure
impulse(LateralSSFull);
title('Full sys Impulse Response')
grid on
% Full Initial
figure
initialValue = eye(5);
dumbSys = ss(ALateralFull,[],CLateralFull,[]);
dumbSys.OutputName = {'v';'p';'r';'phi';'psi'};
for i = 1:5
    initial(dumbSys,initialValue(:,i))
    hold on
end
legend('v','p','r','phi','psi')
title('Full Sys Unity Initial Value Response')
grid on


%% Simulating The Dutch Mode
% Step
figure
step(LateralSSDutch);
title('Dutsh Mode Step Response')
grid on
% Impulse
figure
impulse(LateralSSDutch);
title('Dutsh Mode Impulse Response')
grid on
% Initial Value
figure
initialValue = eye(3);
dumbSys2 = ss(ALateralDutch,[],CLateralDutch,[]);
dumbSys2.OutputName = {'v';'p';'r'};
for i = 1:3
    initial(dumbSys2,initialValue(:,i))
    hold on
end
legend('v','r')
title('Dutsh Mode Unity Initial Value Response')
grid on

%% Simulating The Roling Mode
% Step
figure
step(LateralSSRolling);
title('Rolling Mode Step Response')
grid on
% Impulse
figure
impulse(LateralSSRolling);
title('Rolling Mode Impulse Response')
grid on
% Initial Value
figure
initialValue = eye(1);
dumbSys3 = ss(ALateralRolling,[],CLateralRolling,[]);
dumbSys3.OutputName = {'p'};
initial(dumbSys3,initialValue)
title('Rolling Mode Unity Initial Value Response')
grid on

%% Simulating The Spiral Mode
% Step
figure
step(LateralSSSpiral);
title('Spiral Mode Step Response')
grid on
% Impulse
figure
impulse(LateralSSSpiral);
title('Spiral Mode Impulse Response')
grid on
% Initial Value
figure
initialValue = eye(3);
dumbSys4 = ss(ALateralSpiral,[],CLateralSpiral,[]);
dumbSys4.OutputName = {'p';'r';'phi'};
for i = 1:3
    initial(dumbSys4,initialValue(:,i))
    hold on
end
legend('p','r','phi')
title('Spiral Mode Unity Initial Value Response')
grid on


