%% WRLS estimation of the system parameters and its' modifications
% with variable forgetting factor
% 
% Project prepared for System Identification course
% Author: Bartlomiej Buczek
% 2016
%
% The modifications are covered in related papers:
% http://ieeexplore.ieee.org/document/104093/
% http://ieeexplore.ieee.org/document/6485421/
%%
clear all; close all, clc %#ok<CLSCR,DUALC>

%% Simulation parameters
iter = 700;                                       % iteration number
t = 1:iter;                                       % wektor iteracji

%% System parameters 
% system is defined as y(t) = a(t) * x + b(t)
% a and b are estimated parameters

% a, b definitons
a = [zeros(1,50), linspace(0,1,50), ones(1,100), linspace(1,0,250),...
ones(1,100), zeros(1,100), ones(1,10)./2,zeros(1,40)];

b = [zeros(1,50) ones(1,200) zeros(1,300) ones(1,50) zeros(1,100)];

%% simple WRLS parameters
lambda_1 = 0.9;

%% first WRLS modification parameters
lambda_2 = 1;                                     % lambda zero condition
lambda2_min = 0;

%% second WRLS modification parametrs
% according to paper
lambda_3 = 1;                                     % lambda zero condition
L = 2;
Ka = 2;
Kb = 18;
alpha = 1 - (1 / (Ka * L));
beta = 1 - (1 / (Kb * L));
gamma = 1.2;
psi = 10e-8;
% other zero conditions
vare = 0;
varv = 0;
qkr = 10;

%% Test signal and system response generation
std_in = 1;                                       % test sig std deviation
sig_test = std_in * randn(1,iter);                % test signal (random)

std_mea = 1e-4;                                   % measurment sys std dev.
sig_out_ideal = sig_test .* a + b;                % ideal system response
sig_out_real = sig_test .* a...                   % real system response
    + b + std_mea * randn(1,iter);

%% Matrices initialization
A_1=[]; A_2=[]; A_3=[];
B_1=[]; B_2=[]; B_3=[];
LAMBDA_2=[]; LAMBDA_3=[];
SIG_OUT_1=[]; SIG_OUT_2=[]; SIG_OUT_3=[];

%% Simple WRLS estimation
% zero conditions
est = [0; 0];
P = eye(2) .* 10e3;

for n = 1:iter
    c = ([sig_test(n); 1]' * P * [sig_test(n); 1] + lambda_1)^-1;
    K = P * [sig_test(n); 1] * c;
    P = 1 / lambda_1 .* (P - K * [sig_test(n); 1]' * P);
    % estimation error
    e = sig_out_real(n) - [sig_test(n); 1]' * est;
    est = est + K * e; 
    % store results
    A_1 = [A_1 est(1)];
    B_1 = [B_1 est(2)];
    SIG_OUT_1 =  [SIG_OUT_1 e];
end

%% 1st WRLS modification
% zero conditions
est = [0; 0];
P = eye(2) .* 10e3;

for n = 1:iter
    alfa = sig_out_real(n) - [sig_test(n); 1]' * est;
    L = -round((20*(alfa).^2));
    % Calculate forgetting factor for n-th sample
    lambda_2 = lambda2_min + (1 - lambda2_min) .* (2 .^ L);
    lambda_2 = min(lambda_2, 1);
    
    c = ([sig_test(n); 1]' * P * [sig_test(n); 1] + lambda_2)^-1;
    K = P * [sig_test(n); 1] * c;          
    P = 1 / lambda_2 .* (P - K * [sig_test(n); 1]' * P);
    % Estimation error
    e = sig_out_real(n) - [sig_test(n); 1]' * est;
    est = est + K * e;
    % store results
    A_2 = [A_2 est(1)];
    B_2 = [B_2 est(2)];
    LAMBDA_2 = [LAMBDA_2 lambda_2];
    SIG_OUT_2 =  [SIG_OUT_2 e];
end

%% 2nd WRLS odification
% zero conditions
est = [0; 0];
P = eye(2) .* 10e3;

for n = 1:iter
    q = [sig_test(n); 1]' * P * [sig_test(n); 1];
    qkr = alpha * qkr + (1 - alpha) * q;
    e = sig_out_real(n) - [sig_test(n); 1]' * est;
    vare = alpha * vare + (1 - alpha) * e ^ 2;
    varv = beta * varv + (1 - beta) * e ^ 2;
    % lambda calculation
    if sqrt(vare) <= (gamma * sqrt(varv));
        lambda_3 = 0.999;
    else
        lambda_3 = min(qkr * varv / (psi + abs(vare - varv)), 0.999);
    end

    c = ([sig_test(n); 1]' * P * [sig_test(n); 1] + lambda_3)^-1;
    K = P * [sig_test(n); 1] * c;
    P = 1 / lambda_3 .* (P - K * [sig_test(n); 1]' * P);
    % estimation error
    e = sig_out_real(n) - [sig_test(n); 1]' * est;
    est = est + K * e;
    % Store results
    A_3 = [A_3 est(1)];
    B_3 = [B_3 est(2)];
    LAMBDA_3 = [LAMBDA_3 lambda_3];
    SIG_OUT_3 =  [SIG_OUT_3 e];
end

%% Errors calculation
% simple error
delta_a1 = a - A_1;    delta_a2 = a - A_2;    delta_a3 = a - A_3;
delta_b1 = b - B_1;    delta_b2 = b - B_2;    delta_b3 = b - B_3;

%
D2_a1 = delta_a1.^2;   D2_a2 = delta_a2.^2;   D2_a3 = delta_a3.^2; 
D2_b1 = delta_b1.^2;   D2_b2 = delta_b2.^2;   D2_b3 = delta_b3.^2;

% Mean squared error
error_a1= mean(D2_a1);      error_a2= mean(D2_a2);      error_a3= mean(D2_a3); 
error_b1= mean(D2_b1);      error_b2= mean(D2_b2);      error_b3= mean(D2_b3);

% sum errors of a and b parameters estimation
error_WRLS = error_a1 + error_b1;
error_WRLS_mod1 = error_a2 + error_b2;
error_WRLS_mod2 = error_a3 + error_b3;

%% Results visualization
figure(1)
subplot(211);
plot(t, a,'r',t, A_1,'k',t, A_2,'b',t, A_3,'g',t, a,'r');
title('Parameter a estimation'); grid on;
legend('a real value','WRLS','No 1 modification ', 'No 2 modification');
ylabel('Amplitude'); xlabel('samples');

subplot(212);
plot(t, b,'r',t, B_1,'k',t, B_2,'b',t, B_3,'g',t, b,'r');
title('Parameter b estimation'); grid on;
legend('b real value','WRLS','No 1 modification ', 'No 2 modification');
ylabel('Amplitude'); xlabel('samples');

figure(2)

subplot(211);
plot(t, LAMBDA_2); grid on;
legend('lambda'); grid on; ylabel('lambda'); xlabel('samples');
title('Change in forgetting factor - 1st modification');

subplot(212);
plot(t, LAMBDA_3); grid on;
legend('lambda'); grid on; ylabel('lambda'); xlabel('samples')
title('Change in forgetting factor - 2nd modification')


