%% Steinmetz Constant slip
%Matthew Shirley
%Wednesday evening
%This is consistant with the writeup on overleaf

%Solves for current in stator I_s, given constant slip s.
%Produces a plot of current vs. time and a spectrum, both dimensional

%HEALTH WARNING!!!!! SPECTRUM WILL BE WRONG, need to talk to someone about
%details of evaluating it e.t.c..
% <--- I think this is done correctly now, MRM. 22/07/20

clear all, close all, clc

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

%% calculate dimensional parameter groupings

%The dimensional coefficients:
%a_2_dim \ddot{I}_s + a_1_dim \dot{I}_s 
%                         + a_0_dim I_s =c_1_dim\dot{V}_s +c_0_dim V_s

%These are quantities \tilde{a}_i and \tilde{c}_i on overleaf
a_2_dim = X_r_dash + X_r_dash * X_s / X_m + X_s;

a_1_dim = X_r_dash * R_s/X_m + R_r_dash * X_s/(s*X_m) + R_s + R_r_dash/s;

a_0_dim = R_r_dash * R_s /(s*X_m);

c_1_dim = X_r_dash/X_m + 1;

c_0_dim = R_r_dash/(s*X_m);

%% Calculate NonDim coeffefiecents and define NonDim voltage functions

%Define voltage as funciton handle, varying with time
V_s = @(t) sqrt(2) * cos(2*pi*t); %The nondimensional voltage V_s(t)

dV_s = @(t) - 2 * sqrt(2) * pi * sin(2*pi*t); %The nondimensional derivative of voltage \dot{V}_s(t)

%Nondimensional coefficients:
% \ddot{I}_s + a_1 \dot{I}_s + a_2 \dot{I}_s = c_1 \dot{V}_s + c_0 V_s

a_1 = (a_1_dim/a_2_dim) * t0;

a_0 = (a_0_dim/a_2_dim)  * t0^2;

c_1 = (c_1_dim/a_2_dim) * (V0/I0) * t0; 

c_0 = (c_0_dim/a_2_dim) * (V0/I0) * t0^2; 

%% Define First order ODE system and solve for I_s
%Solve ODE using ode45

%ode45 only accepts first order system. 
%Let u = [I_s, dI_s]

%Define RHS of first order system
RHS = @(t,u) [0, 1; - a_0, -a_1]* u  + [0; c_1 * dV_s(t) + c_0 *V_s(t)];

%Define options for ODE solver
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

%Define end time
t_end = 50;
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