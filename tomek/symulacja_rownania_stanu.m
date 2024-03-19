clear all;

run("stale.m");

simulation_time = 2500; % Czas symulacji (sekundy)
time_step = 0.1; % Krok symulacji (sekundy)

tau_steps = tau / time_step;
tau_C_steps = tau_C / time_step;
k_min = max(tau_steps, tau_C_steps) + 1;
k_max = simulation_time / time_step + k_min;

% Trajektorie wejść procesu
T_H(1:k_max) = T_Hpp;
T_C(1:k_max) = T_Cpp;
T_D(1:k_max) = T_Dpp;
F_H(1:k_max) = F_Hpp;
F_Cin(1:k_max) = F_Cpp;
F_D(1:k_max) = F_Dpp;

F_H(500/time_step:k_max) = F_Hpp * 1.1;
F_H(1000/time_step:k_max) = F_Hpp * 0.9;
F_Cin(1500/time_step:k_max) = F_Cpp * 1.1;
F_Cin(2000/time_step:k_max) = F_Cpp * 0.9;

% Stan i wyjścia procesu przed rozpoczęciem symulacji
x_1(1:k_min) = C * h_pp^3;
x_2(1:k_min) = C * h_pp^3 * T_pp;
h(1:k_min) = h_pp;
T(1:k_min) = T_pp;
T_out(1:k_min) = T_pp;

% Symulacja (na podstawie równań stanu)
for k = k_min:k_max
    F_C = F_Cin(k-tau_C_steps);
    u_1 = F_H(k);
    u_2 = F_C;
    u_3 = F_D(k);
    u_4 = F_H(k) * T_H(k);
    u_5 = F_C * T_C(k);
    u_6 = F_D(k) * T_D(k);

    dx_1 = -0.0068909048699667508129585103576364 * x_1(k-1) + u_1 + u_2 + u_3 - 44.999999999999999953822964851574;
    dx_2 = 1.2689601318043770646157656936667 * x_1(k-1) - 0.041345429219800504877751062145819 * x_2(k-1) + u_4 + u_5 + u_6 - 1657.349999999999870839276406354;
    y_1 = 0.0033076343375840403868259350368931 * x_1(k-1) + 8.6399999999999999822680185030046;
    y_2 = -0.028199114040097268131509270530017 * x_1(k-1) + 0.00076565609666297231333662845056432 * x_2(k-1) + 36.829999999999997167555035797123;
    x_1(k) = x_1(k-1) + time_step * dx_1;
    x_2(k) = x_2(k-1) + time_step * dx_2;

    h(k) = y_1;
    T(k) = y_2;
    T_out(k) = T(k-tau_steps);
end

figure(2);
plot(T_out(k_min:k_max));
hold on
plot(h(k_min:k_max));
hold off
legend(["T_{out}", "h"]);



