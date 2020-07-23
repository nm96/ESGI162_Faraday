%% Varaible slip parameters
%Matthew Shirley
%Wednesday

%% Define dinmensional parameters

freq = 50;  %frequency Hertz

s = 0.05; %Slip 

R_s = 2.3; %stator resistance (Ohms)

X_s =  17.5e-3 ;  %Stator reactance (Ohms)

R_r_dash = 4.2; %Rotor resistance (Ohms)

X_r_dash = X_s;  %Leakage resistance (Ohms)


X_m = 72.25/(2 * pi * freq);  %Magnetising resistance (Ohms)

V_rms = 90; %root mean squared of voltage


%% Pick nondimensional scalings
%Pick the scalings used to nondimensionalise the model

t0 = 1/freq;   %timescale  t= t0 \hat{t}
V0 = V_rms;    %voltage    V= V_{rms} \hat{V}
I0 = 1/R_r_dash; %current    I = \frac{V_{rms}}{R_r_dash}

%% Calcualate dimensional parameter groupings 

%Write dimensional equation as:
% ( \tilde{a}_{2,1}s ) \ddot{I}_s 
% + ( \tilde{a}_{1,2}\dot{s} + \tilde{a}_{1,1}s + \tilde{a}_{1,0} ) \dot{I}_s 
% + ( \tilde{a}_{0,2}\dot{s} + \tilde{a}_{0,0} ) I_s 
%= ( \tilde{c}_{1,1}s )\dot{V}_s + ( \tilde{c}_{0,2}\dot{s} + \tilde{c}_{0,0} )V_s.


%Note that a_2_0_dim = \tilde{a}_{2,0} in Overleaf.

a_2_1_dim = X_r_dash + X_r_dash * X_s /X_m + X_s;

a_1_2_dim = X_r_dash + X_r_dash * X_s /X_m + X_s;

a_1_1_dim = X_r_dash * R_s /X_m + R_s;

a_1_0_dim = R_r_dash + R_r_dash * X_s / X_m;

a_0_2_dim = X_r_dash * R_s /X_m + R_s;

a_0_0_dim = R_r_dash * R_s / X_m;

c_1_1_dim = X_r_dash/X_m + 1;

c_0_2_dim = X_r_dash/X_m + 1;

c_0_0_dim = R_r_dash/X_m;

%% Calculate nondimensional parameter groupings

%Write nondimensional equation as:
% s \ddot{I}_s
%  + ( a_{1,2}\dot{s} + a_{1,1}s + a_{1,0} ) \dot{I}_s 
%        + ( a_{0,2}\dot{s} + a_{0,0} ) I_s 
%           =  ( c_{1,1}s)\dot{V}_s  + ( c_{0,2}\dot{s} + c_{0,0} ),

a_1_2 = a_1_2_dim/ a_2_1_dim

a_1_1 = a_1_1_dim / a_2_1_dim * t0

a_1_0 = a_1_0_dim / a_2_1_dim * t0

a_0_2 = a_0_2_dim / a_2_1_dim * t0

a_0_0 = a_0_0_dim / a_2_1_dim * t0^2

c_1_1 = c_1_1_dim / a_2_1_dim  * V0/I0 * t0

c_0_2 = c_0_2_dim/ a_2_1_dim * V0/I0 * t0

c_0_0 = c_0_0_dim/ a_2_1_dim * V0/I0 * t0^2

%Define voltage as funciton handle, varying with time
V_s = @(t) sqrt(2) * cos(2*pi*t); %The nondimensional voltage V_s(t)

dV_s = @(t) - 2 * sqrt(2) * pi * sin(2*pi*t); %The nondimensional derivative of voltage \dot{V}_s(t)


% define slip and its derivative
omega = 2*pi;

s= @(t)  0.05 + 0.005 * cos(omega * t);
%s = @(t) 0.05;

ds = @(t)  - 0.005 * omega * sin(omega * t);
%ds = @(t) 0;

%% Rewrite as a first order system

RHS_2_1 =  @(t) - ( a_0_2 * ds(t)/s(t)  + a_0_0/s(t));

RHS_2_2 = @(t)  - ( a_1_2  * ds(t)/s(t)  + a_1_1  + a_1_0 /s(t));

RHS = @(t, u) [0, 1; RHS_2_1(t), RHS_2_2(t)]* u ...
    + [0; c_1_1 * dV_s(t) + (c_0_2 * ds(t)/s(t) + c_0_0 /s(t)) *V_s(t)];

%Define options for ODE solver
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

%Define end time
t_end = 2000;
tspan = [0, t_end]; %Interval to solve over

%inital condition
u_inital = [0; 0];

%solve using ode45
[t, u] = ode45(RHS, tspan, u_inital,options);

%t gives vector of timesteps
% u = [I_s, dI_s/dt] evaluated at points given by t

%% Convert solution to dimensional variables

%Dimensional time
t_dim = t0 * t;
%Dimensional current I_s
I_s_dim = I0 * u(:,1);

%% Plot current vs. time
ax1 = axes;
plot(ax1,t_dim,I_s_dim,'b-','LineWidth',1.0)
ax1.FontSize = 14;
xlabel('time, $\mathrm{sec}$','Interpreter','latex')
ylabel('Current in stator, $\mathrm{A}$','Interpreter','latex')


%% Plot spectrum

% Calculate the fast Fourier Transform  of I_s
% MRM - I think this does the fft bit for us.

Isp = spline(t_dim,I_s_dim);

SampleFreq = 1000;
SamplePeriod = 1/SampleFreq;
SigLength = max(t_dim)/SamplePeriod;

T = (0:1:SigLength-1)*SamplePeriod;

Ispdata = ppval(Isp,T);
Ispdata(1) = Ispdata(1)/2 + Ispdata(end)/2;

spec = fft(Ispdata);
spec = abs(spec(1:SigLength/2+1)/SigLength); % Frequencies less than Nyquist

f = SampleFreq*(0:(SigLength/2))/SigLength;

%Define range
freq_min = 0;
freq_max = 100;

%plot
figure
ax2 =axes;
plot(f,spec,'b-o','LineWidth',0.8)
xlim([freq_min, freq_max])
xlabel('frequency, $\mathrm{Hz}$','Interpreter','latex')