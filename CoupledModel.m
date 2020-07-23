%% Script to solve the combined model.
%Matthew Shirley

%Soves for u, v, and theta

%% Define parameters 

M = 1

c_x = 0.1;
c_y = 0.1;
c_theta = 0.1;

k_x = 100;
k_y = 100;

epsilon = 0.01;

I = 1

Phi = 1
R = 1
N = 24
Omega = 50

contact_1 = @(t,w) 0;
contact_2 = @(t,w) 0;
contact_3 = @(t,w) 0;


%% Define components of first order system 
%Rewrite equations as first order system (see overleaf)

% Will solve for  = [u, \dot{u}, v, \dot{v}, \theta, \dot{\theta}]

%Define 6 * 6 mass matrix
Mass = @(t,w) [ 0, M, 0, 0, 0, -epsilon * M * sin(w(5)) ; ...
            0, 0, 0, M, 0, epsilon * M * cos(w(5)); ...
            0, M * w(3) - epsilon * M * sin(w(5)) , 0 , - M*w(1) + epsilon * M * cos(w(5)), (I + M*epsilon^2) - epsilon * M * w(1) * cos(w(5)) - epsilon * M * w(3) * sin(w(5)) ; ...
            1, 0, 0, 0, 0, 0; ...
            0, 0, 1, 0, 0, 0; ...
            0, 0, 0, 0, 1, 0];


% RHS 

RHS_1 = @(t,w)  epsilon * M * (w(6))^2 * cos(w(5)) - c_x *w(2) - k_x *w(1) + contact_1(t,w);

RHS_2 = @(t,w) epsilon * M * (w(6))^2 * sin(w(5)) - c_y * w(4) - k_y * w(3) + contact_2(t,w);

RHS_3 = @(t,w) - Phi^2 / R * N/2 * w(5) +  Phi^2 / R * N/2 * Omega  -c_theta * w(6) ...
    + contact_3(t,w) - epsilon * M * w(1) * (w(6))^2 * sin(w(5)) + epsilon * M * w(3) * (w(6))^2 * cos(w(5));

RHS_4 = @(t,w) w(2);

RHS_5 = @(t,w) w(4);

RHS_6 = @(t,w) w(6);

RHS = @(t,w) [ RHS_1(t,w) ; RHS_2(t,w); RHS_3(t,w); RHS_4(t,w); RHS_5(t,w); RHS_6(t,w)];


%% Define ODE solver

%Define options for ODE solver
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Mass',M);

%Define end time
t_end = 100;
tspan = [0, t_end]; %Interval to solve over

%initial condition
w0 = [0; 0; 0; 0; 0; Omega];

%solve using ode45
[t, w] = ode45(RHS, tspan, w0,options);

%extract components
u = w(:,1);
v = w(:,3);
theta = w(:,5);
dtheta = w(:,6);
%% plot results

figure
ax1 = axes;
plot(ax1,t,u,'b-','LineWidth',1.0)
ax1.FontSize = 14;
xlabel('time','Interpreter','latex')
ylabel('u','Interpreter','latex')

figure
ax2 = axes;
plot(ax2,t,v,'b-','LineWidth',1.0)
ax2.FontSize = 14;
xlabel('time','Interpreter','latex')
ylabel('v','Interpreter','latex')

figure
ax3 = axes;
plot(ax3,t,theta,'b-','LineWidth',1.0)
ax3.FontSize = 14;
xlabel('time','Interpreter','latex')
ylabel('theta','Interpreter','latex')

%% Calculate current

%HEALTH WARNING !!!!!!!!
%Use equation at end of section 2.3 on overleaf to calculate I_{res}
% Not sure this equation consistant with ones I used because has \dot{Phi}
% in it and some subscript i 

A_s = 1;
L_s = 1;
k=1;
phi_s = 1;
dPhi = 1;

I_res = A_s /(2 *L_s) * k/R * ( Phi *(Omega - dtheta).*sin(Omega * t - phi_s) + dPhi * cos(Omega * t - phi_s));

%% Obtain spectrum of I_res

%HEALTH WARNING may need to redimensionalise
%Copied from Matt M Steintmetz code

% Calculate the fast Fourier Transform  of I_s
% MRM - I think this does the fft bit for us.

Isp = spline(t,I_res);

SampleFreq = 1000;
SamplePeriod = 1/SampleFreq;
SigLength = max(t)/SamplePeriod;

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
semilogy(f,spec,'b-o','LineWidth',0.8)
xlim([freq_min, freq_max])
xlabel('frequency, $\mathrm{Hz}$','Interpreter','latex')