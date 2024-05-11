clear all;
close all;
% clc;

% vector = [-1.25, -0.005, -0.04, 1.25, 0.005, 0.004];
vector = [-1.05, 200, 0.01, 1.05, 100, 0.01];
%vector = [-1.05, 200, 0.01, 1.05, 100, 0.01];
Kp_1=vector(1); % człon proporcjonalny
Ki_1=vector(2); % człon całkujący
Kd_1=vector(3);
Kp_2=vector(4); % człon proporcjonalny
Ki_2=vector(5); % człon całkujący
Kd_2=vector(6);

run("stale.m");

simulation_time = 10000; % Czas symulacji (sekundy)
T_p = 1; % Krok symulacji (sekundy)

tau_steps = tau / T_p;
tau_C_steps = tau_C / T_p;
k_min = max(tau_steps, tau_C_steps) + 1;
k_max = simulation_time / T_p + k_min;

% Trajektorie wejść procesu
T_H(1:k_max) = T_Hpp;
T_C(1:k_max) = T_Cpp;
T_D(1:k_max) = T_Dpp;
F_H(1:k_max) = F_Hpp;
F_Cin(1:k_max) = F_Cpp;
F_D(1:k_max) = F_Dpp;

u1 = F_Cin;
u2 = F_H;

T_zad(1:k_max) = 36.83;
T_zad(round(k_max/3):k_max) = 36.83 - 2;
T_zad(round(2*k_max/3):k_max) = 36.83 + 2;

h_zad(1:k_max) = 12.96;
h_zad(round(k_max/3):k_max) = 12.96 + 1;
h_zad(round(2*k_max/3):k_max) = 12.96 - 1;

% Stan i wyjścia procesu przed rozpoczęciem symulacji
F(1:k_max) = alpha * sqrt(h_pp);
F_C(1:k_max) = F_Cpp;
h(1:k_max) = h_pp;
V(1:k_max) = C * h_pp^3;
T(1:k_max) = T_pp;
VT(1:k_max) = V(1) * T(1);
T_out(1:k_max) = T_pp;

e1(1:k_max) = 0;
e2(1:k_max) = 0;
e = 0;

r0_y1 = Kp_1*(T_p/(2*Ki_1) + Kd_1/T_p + 1);
r1_y1 = Kp_1*(T_p/(2*Ki_1) - 2*Kd_1/T_p - 1);
r2_y1 = Kp_1*Kd_1/T_p;

r0_y2 = Kp_2*(T_p/(2*Ki_2) + Kd_2/T_p  + 1);
r1_y2 = Kp_2*(T_p/(2*Ki_2) - 2*Kd_2/T_p - 1);
r2_y2 = Kp_2*Kd_2/T_p;

for k = k_min:k_max
 %% bez odsprzęgania
 %  [F_C, V, VT, T, F, h, T_out] = obiekt(F_Cin, F_H, F_D, F_C, T_H, T_C, T_D, T_out, h, C, alpha, tau_C_steps, tau_steps, V, VT, T, F, T_p, k);
 %%
    e1(k) = T_zad(k) - T_out(k-1);
    F_Cin(k) = F_Cin(k-1) + r0_y1*e1(k) + r1_y1*e1(k-1) + r2_y1*e1(k-2);
    if F_Cin(k) < 0
        F_Cin(k) = 0;
    end
    
    e2(k) = h_zad(k) - h(k-1);
    F_H(k) = F_H(k-1) + r0_y2*e2(k) + r1_y2*e2(k-1) + r2_y2*e2(k-2);
    if F_H(k) < 0
        F_H(k) = 0;
    end

    e = e + abs(e1(k))^2 + abs(e2(k))^2;

   
     %% z odsprzęganiem
    u1(k) = F_Cin(k-1) - 0.2*(F_H(k-1) - F_Hpp);
    u2(k) = F_H(k-1) - 0.3*(F_Cin(k-1) - F_Cpp);

    if u1(k) < 0
        u1(k) = 0;
    end
    if u2(k) < 0
        u2(k) = 0;
    end
    [F_C, V, VT, T, F, h, T_out] = obiekt(u1, u2, F_D, F_C, T_H, T_C, T_D, T_out, h, C, alpha, tau_C_steps, tau_steps, V, VT, T, F, T_p, k);
    
end
disp(e)

figure(1)
hold on
plot(T_out(k_min:k_max))
plot(T_zad(k_min:k_max))
plot(h(k_min:k_max))
plot(h_zad(k_min:k_max))
hold off
legend(["T_{out}", "Tzad", "h", "hzad"], Location="best")

figure(2)
hold on
plot(F_Cin(k_min:k_max))
plot(F_H(k_min:k_max))
hold off
legend(["Fcin", "Fh"])
