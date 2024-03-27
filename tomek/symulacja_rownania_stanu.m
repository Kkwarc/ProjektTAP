clear all;

run("stale.m");

simulation_time = 2500; % Czas symulacji (sekundy)
T_p = 0.1; % Krok symulacji (sekundy)

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

F_H(500/T_p:k_max) = F_Hpp * 1.1;
F_H(1000/T_p:k_max) = F_Hpp * 0.9;
F_Cin(1500/T_p:k_max) = F_Cpp * 1.1;
F_Cin(2000/T_p:k_max) = F_Cpp * 0.9;

% Stan i wyjścia procesu przed rozpoczęciem symulacji
x_1 = C * h_pp^3;
x_2 = C * h_pp^3 * T_pp;
x = [x_1, x_2]';
h(1:k_min) = h_pp;
T(1:k_min) = T_pp;
T_out(1:k_min) = T_pp;

% Wczytanie macierzy A, B, C, D
run("wspolczynniki_rownania_stanu.m");

% Symulacja (na podstawie równań stanu)
% Metoda RK2
for k = k_min:k_max
    F_C = F_Cin(k-tau_C_steps);
    u_1 = F_H(k);
    u_2 = F_C;
    u_3 = F_D(k);
    u_4 = F_H(k) * T_H(k);
    u_5 = F_C * T_C(k);
    u_6 = F_D(k) * T_D(k);
    u = [u_1, u_2, u_3, u_4, u_5, u_6, 1]';

    k_1 = A*x + B*u;
    x_temp = x + 0.5 * T_p * k_1;
    k_2 = A*x_temp + B*u;
    x = x + T_p * k_2;
    y = C*x + D*u;

    h(k) = y(1);
    T(k) = y(2);
    T_out(k) = T(k-tau_steps);
end

figure(7);
hold on;
time = (k_min-1:k_max-1) * T_p;
plot(time, T_out(k_min:k_max));
hold on
plot(time, h(k_min:k_max));
hold off
legend(["T_{out}", "h"]);



