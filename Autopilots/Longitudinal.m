clear; clc; close all;

%% Reading Aircraft Data From Excel
Aircraft_data = xlsread('B-747_data.xlsx','B1:B50');

%% Flight Conditoins Data
flight = num2cell(Aircraft_data(1:17));
[W,g,s,b,c,h,M,V_tot,CG,Ix,Iy,Iz,Ixz,Ixy,Iyz,Alpha,q]= deal(flight{:});
m = W/g;
%% Longitudinal Dinemntional Stability Derivatives Data
long = num2cell(Aircraft_data(18:35));
[u,w,v,Xu,Zu,Mu,Xw,Zw,Mw,Xq,Zq,Mq,XDE,ZDE,MDE,XDTH,ZDTH,MDTH] = deal(long{:});

%% Longitudinal non-Dinemntional Stability Derivatives Data
Cxu = Xu*m*u/(q*s);
Czu = Zu*m*u/(q*s);
Cmu = Mu*Iy*u/(q*s*c);
Cxw = Xw*m*u/(q*s);
Czw = Zw*m*u/(q*s);
Cmw = Mw*Iy*u/(q*s*c);
Cxq = 2*Xq*m*u/(q*s*c);
Czq = 2*Zq*m*u/(q*s*c);
Cmq = 2*Mq*Iy*u/(q*s*c^2);
Cxe = XDE*m/(q*s);
Cze = ZDE*m/(q*s);
Cme = MDE*Iy/(q*s*c);
Cxth = XDTH*m/(q*s);
Czth = ZDTH*m/(q*s);
Cmth = MDTH*Iy/(q*s*c);

%% Longitudinal State Space Model
ALongFull = [Xu Xw Xq -g; Zu Zw Zq+u 0; Mu Mw Mq 0; 0 0 1 0];
BLongFull = [XDE XDTH; ZDE ZDTH; MDE MDTH; 0 0];
CLongFull = eye(4);
DLongFull = zeros(4,2);
LongSSFull = ss(ALongFull,BLongFull,CLongFull,DLongFull);
LongSSFull.OutputName = {'u';'w';'q';'theta'};
LongSSFull.InputName = {'De';'Dth'};

%% Longitudinal Short & Long Period Modes
% Short Mode
ALongShort = [Zw Zq+u;Mw Mq];
BLongShort = [ZDE ZDTH; MDE MDTH];
CLongShort = eye(2);
DLongShort = zeros(2,2);
LongSSShort = ss(ALongShort,BLongShort,CLongShort,DLongShort);
LongSSShort.OutputName = {'w';'q'};
LongSSShort.InputName = {'De';'Dth'};
% Long Mode
ALongLong = [Xu -g;-Zu/(u+Zq) 0];
BLongLong = [XDE XDTH; -ZDE/(u+Zq) -ZDTH/(u+Zq)];
CLongLong = eye(2);
DLongLong = zeros(2,2);
LongSSLong = ss(ALongLong,BLongLong,CLongLong,DLongLong);
LongSSLong.OutputName = {'u';'theta'};
LongSSLong.InputName = {'De';'Dth'};

%% Pole Placement In Case Of Insatiability
PolesFull = pole(LongSSFull);
if (isreal(sqrt(real(PolesFull))))  % enters iff there is a positive pole
    % Full
    PolesFullNew = complex(-abs(real(PolesFull)),imag(PolesFull));
    kFull = place(ALongFull,BLongFull,PolesFullNew);
    ALongFull = ALongFull-BLongFull*kFull;
    LongSSFull = ss(ALongFull,BLongFull,CLongFull,DLongFull);
    LongSSFull.OutputName = {'u';'w';'q';'theta'};
    LongSSFull.InputName = {'De';'Dth'};
    % Short Mode
    PolesShort = pole(LongSSShort);
    PolesShortNew = complex(-abs(real(PolesShort)),imag(PolesShort));
    kShort = place(ALongShort,BLongShort,PolesShortNew);
    ALongShort = ALongShort-BLongShort*kShort;
    LongSSShort = ss(ALongShort,BLongShort,CLongShort,DLongShort);
    LongSSShort.OutputName = {'w';'q'};
    LongSSShort.InputName = {'De';'Dth'};
    % Long Mode
    PolesLong = pole(LongSSLong);
    PolesLongNew = complex(-abs(real(PolesLong)),imag(PolesLong));
    kLong = place(ALongLong,BLongLong,PolesLongNew);
    ALongLong = ALongLong-BLongLong*kLong;
    LongSSLong = ss(ALongLong,BLongLong,CLongLong,DLongLong);
    LongSSLong.OutputName = {'u';'theta'};
    LongSSLong.InputName = {'De';'Dth'};
end

%% Longitudinal Eigenvalues, eigenvectors
[eigVectorsLongFull, eigValuesLongFull] = eig(ALongFull);

%% Full & Modes Characteristics (undamped natural frequency, damping ratio, settling time, etc.)
% Full Mode
[WnFull,ZetaFull] = damp(LongSSFull);
StepInfoFull = stepinfo(step(LongSSFull));
% Short Mode
[WnShort,ZetaShort] = damp(LongSSShort);
StepInfoShort = stepinfo(step(LongSSShort));
% Long Mode
[WnLong,ZetaLong] = damp(LongSSLong);
StepInfoLong = stepinfo(step(LongSSLong));

%% Modes Approx. Error
ErorrShort = abs(WnShort-WnFull(3:4))/WnFull(3:4)*100;
ErorrLong = abs(WnLong-WnFull(1:2))/WnFull(1:2)*100;

%% Longitudinal Full Transfer Function
LongtfFull = tf(LongSSFull);
LongtfFull_u_de = LongtfFull(1,1);
LongtfFull_u_dT = LongtfFull(1,2);
LongtfFull_w_de = LongtfFull(2,1);
LongtfFull_w_dT = LongtfFull(2,2);
LongtfFull_q_de = LongtfFull(3,1);
LongtfFull_q_dT = LongtfFull(3,2);
LongtfFull_theta_de = LongtfFull(4,1);
LongtfFull_theta_dT = LongtfFull(4,2);

%% Longitudinal Modes Transfer Function
% Short Period Mode
LongtfShort = tf(LongSSShort);
LongtfShort_w_de = LongtfShort(1,1);
LongtfShort_w_dT = LongtfShort(1,2);
LongtfShort_q_de = LongtfShort(2,1);
LongtfShort_q_dT = LongtfShort(2,2);
% Long Period Mode
LongtfLong = tf(LongSSLong);
LongtfLong_u_de = LongtfLong(1,1);
LongtfLong_u_dT = LongtfLong(1,2);
LongtfLong_theta_de = LongtfLong(2,1);
LongtfLong_theta_dT = LongtfLong(2,2);

%% Long Autopilot parameters
servo=tf(10,[1,10]);
engine_timeLag=tf(0.1,[1,0.1]);
HPC=tf([1 0],[1 0.4]);
u_0=u;
W_0=u_0*tan(Alpha);
thetha_0=Alpha;
height_o=0;

%% Simulating The Full System
% Full Step
step(LongSSFull);
title('Full sys Step Response')
grid on
% Full Impulse
figure
impulse(LongSSFull);
title('Full sys Impulse Response')
grid on
% Full Initial
figure
initialValue = eye(4);
dumbSys = ss(ALongFull,[],CLongFull,[]);
dumbSys.OutputName = {'u';'w';'q';'theta'};
for i = 1:4
    initial(dumbSys,initialValue(:,i))
    hold on
end
legend('u','w','q','theta')
title('Full Sys Unity Initial Value Response')
grid on

%% Simulating The Short Mode
% Step
figure
step(LongSSShort);
title('Short Mode Step Response')
grid on
% Impulse
figure
impulse(LongSSShort);
title('Short Mode Impulse Response')
grid on
% Initial Value
figure
initialValue = eye(2);
dumbSys2 = ss(ALongShort,[],CLongShort,[]);
dumbSys2.OutputName = {'w';'q'};
for i = 1:2
    initial(dumbSys2,initialValue(:,i))
    hold on
end
legend('w','q')
title('Short Mode Unity Initial Value Response')
grid on

%% Simulating The Long Mode
% Step
figure
step(LongSSLong,200);
title('Long Mode Step Response')
grid on
% Impulse
figure
impulse(LongSSLong,200);
title('Long Mode Impulse Response')
grid on
% Initial Value
figure
initialValue = eye(2);
dumbSys3 = ss(ALongLong,[],CLongLong,[]);
dumbSys3.OutputName = {'u';'theta'};
for i = 1:2
    initial(dumbSys3,initialValue(:,i),200)
    hold on
end
legend('u','theta')
title('Long Mode Unity Initial Value Response')
grid on


